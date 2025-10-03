
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


#include "eval.h"
#include "param.h"
#include "luna.h"

extern logger_t logger;

extern writer_t writer;

extern freezer_t freezer;


//
// cmd_t
//

cmd_t::cmd_t( const bool silent ) 
{
  register_specials();
  reset();
  error = ! read( NULL , silent );
}

cmd_t::cmd_t( const std::string & str , const bool silent ) 
{
  register_specials();
  reset();
  error = ! read( &str , silent ); 
}

void cmd_t::add_cmdline_cmd( const std::string & c ) 
{
  cmdline_cmds.append( c + " " );
}

void cmd_t::reset() 
{
  cmds.clear();
  params.clear();
  line = "";
  error = false;
  will_quit = false;
}

void cmd_t::clear_static_members() 
{
  input = "";
  cmdline_cmds = "";
  stout_file = "";
  append_stout_file = false;
  
  vars.clear();
  ivars.clear();
  idmapper.clear();
  signallist.clear();
  label_aliases.clear();
  primary_alias.clear();
  primary_upper2orig.clear();
}

bool cmd_t::empty() const 
{ 
  return will_quit; 
}

bool cmd_t::valid() const 
{    
  if ( error ) return false;
  /* for (int c=0;c<cmds.size();c++) */
  /*   if ( commands.find( cmds[c] ) == commands.end() ) return false; */
  return true;
}

bool cmd_t::badline() const 
{ 
  return error; 
} 

std::string cmd_t::offending() const 
{ 
  return ( error ? line : "" ); 
}

int cmd_t::num_cmds() const 
{ 
  return cmds.size(); 
}

std::string cmd_t::cmd(const int i) 
{ 
    return cmds[i]; 
}

param_t & cmd_t::param(const int i) 
{ 
  return params[i]; 
}

bool cmd_t::process_edfs() const
{
  // all commands process EDFs, /except/ the following
  if ( cmds.size() == 1 
       && ( cmds[0] == "" 
	    || cmds[0] == "." 
	    || Helper::iequals( cmds[0] , "DUMMY" ) 
	    || Helper::iequals( cmds[0] , "INTERVALS" )
	    ) )
    return false;
  return true;    
}

bool cmd_t::is( const int n , const std::string & s ) const
{
  if ( n < 0 || n >= cmds.size() ) Helper::halt( "bad command number" );
  return Helper::iequals( cmds[n] , s );
}
  
std::string cmd_t::data() const 
{ 
  return input; 
} 

bool cmd_t::quit() const 
{ 
  return will_quit; 
}

  
void cmd_t::quit(bool b) 
{ 
  will_quit = b; 
}


void cmd_t::signal_alias( const std::string & s )
{

  // the primary alias can occur multiple times, and have multiple 
  // labels that are mapped to it

  // label_aliases[ ALIAS ] -> primary;
  // primary_alias[ primary ] -> [ ALIASES ] 
  // primary_upper2orig[ PRIMARY ] -> primary
  
  // however: two rules
  // 1. many-to-one mapping means the same label cannot have multiple primary aliases
  // 2. following, and on principle of no transitive properties, alias cannot have alias

  // X|Y|Z    
  // X|A|B    okay
  
  // W|A  bad, A already mapped
  // V|X  bad, X already mapped
  // i.e. things can only occur once in the RHS, but multiple times in the LHS

  // keep all RHS aliases as UPPERCASE
  // but with case-insensitive matches
  
  // x|Y|Y

  // format canonical|alias1|alias2 , etc.
  std::vector<std::string> tok = Helper::quoted_parse( s , "|" );    
  if ( tok.size() < 2 ) Helper::halt( "bad format for signal alias:  canonical|alias 1|alias 2\n" + s );
  const std::string primary = Helper::unquote( tok[0] );

  // has the LHS primary already been an alias?
  if ( label_aliases.find( Helper::toupper( primary ) ) != label_aliases.end() )
    Helper::halt( primary + " specified as both primary alias and mapped term" );

  for (int j=1;j<tok.size();j++) 
    {

      // impose rules, use upper case version of all aliases
      const std::string mapped = Helper::unquote( tok[j] );
      const std::string uc_mapped = Helper::toupper( mapped );
      
      if ( primary_upper2orig.find( uc_mapped ) != primary_upper2orig.end() )
	Helper::halt( mapped + " specified as both primary alias and mapped term" );

      // same alias cannot have multiple, different primaries
      // although we track the case of the primary, do case-insensitive match
      
      if ( label_aliases.find( uc_mapped ) != label_aliases.end() )
	if ( ! Helper::iequals( primary , label_aliases[ uc_mapped ] ) )
	  Helper::halt( mapped + " specified twice (case-insensitive) in alias file w/ different primary aliases" );
      
      // otherwise, set this alias, using UC version of the mapped term
      label_aliases[ uc_mapped ] = primary;
      
      primary_alias[ primary ].push_back( uc_mapped );

      // also, for lookup/checks, track UC primary for case-insensitive matches
      // but first check, if we have seen the primary before, it needs to be identical w.r.t. case
      if ( primary_upper2orig.find( Helper::toupper( primary ) ) != primary_upper2orig.end() )
	{
	  if ( primary_upper2orig[ Helper::toupper( primary ) ] != primary )
	    Helper::halt( "primary alias specified with varying case:"
			  + primary_upper2orig[ Helper::toupper( primary ) ]  + " and " + primary );
	    }
      else
	primary_upper2orig[ Helper::toupper( primary ) ] = primary;
    }
    
}


const std::set<std::string> & cmd_t::signals() 
{ 
  return signallist; 
}
  
void cmd_t::clear_signals() 
{ 
  signallist.clear(); 
}

std::string cmd_t::signal_string() 
{
  
  if ( signallist.size() == 0 ) return "*"; // i.e. all signals
  
  std::stringstream ss;
  std::set<std::string>::iterator ii = signallist.begin();
  while ( ii != signallist.end() )
    {
      if ( ii != signallist.begin() ) ss << ",";
      ss << *ii;	  
      ++ii;
    }
  return ss.str();
}
  

void cmd_t::populate_commands() 
{
  // redundant... now using cmddefs_t  
}



// ----------------------------------------------------------------------------------------
//
// Process commands from STDIN
//
// ----------------------------------------------------------------------------------------

void cmd_t::dynamic_parse( const std::string & id )
{
  // new version of cmd_t::replace_wildcards()
  // which dynamically goes through and does
  //   - variable substitutions
  //   - loop unrolling
  //   - block conditionals
  // in a single pass;

  // no commands are evaluated, but this is called per-person to allow generic
  // dynamic specification of scripts: i.e. variables cannot depend on the outputs 
  // of luna commands;

  // dynamic_parse() should be a drop-in replacement for replace_wildcards()
  //  e.g. called from lunapi_t too
  //
  // all the main work is done in proc_all()
  //

  //
  // main input == line (std::string)
  //   - will contain the original script
  //   - comments, newlines versus '&' have been taken care of
  //
  
  
  //
  // If the script contains an ID wildcard, confirm that ID does not already contain the wildcard 
  //
  
  if ( line.find( globals::indiv_wildcard ) != std::string::npos && 
       id.find( globals::indiv_wildcard ) != std::string::npos ) 
    Helper::halt( "ID " + id + " contains ID-wildcard character " 
		  + globals::indiv_wildcard + " (i.e. use wildcard=X to specify in a different one)" );
  

  //
  // Copy a set of variables, where any i-vars will overwrite an existing var
  //
  
  std::map<std::string,std::string> allvars = vars;
  if ( ivars.find( id ) != ivars.end() )
    {
      const std::map<std::string,std::string> & newvars = ivars.find( id )->second;
      std::map<std::string,std::string>::const_iterator vv = newvars.begin();      
      while ( vv != newvars.end() )
	{
	  allvars[ vv->first ] = vv->second;
	  ++vv;
	}
    }  

  
  //
  // Parse into lines
  //

  std::vector<std::string> tok = Helper::quoted_parse( line , "\n" );
    
  
  //
  // expand loops, conditionals and varaibles
  //
  
  tok = cmd_t::proc_all( tok , &allvars );

  //
  // command(s): do this just for printing; real parsing will include variables, etc
  //
  
  params.clear();
  cmds.clear();

  for (int c=0;c<tok.size();c++)
    {      
      std::vector<std::string> ctok = Helper::quoted_parse( tok[c] , "\t " );
      
      // may be 0 if a line was a variable declaration
      if ( ctok.size() >= 1 ) 
	{
	  //std::cout << "adding [" << ctok[0] << "]\n";
	  cmds.push_back( ctok[0] );
	  
	  param_t param;
	  for (int j=1;j<ctok.size();j++) param.parse( ctok[j] );
	  params.push_back( param );
	}
    }
  

  // replace in all 'params' any instances of 'globals::indiv_wildcard' with 'id'
  // ALSO, will expand any @{includes} from files (which may contain ^ ID wildcards
  // in their names,  e..g.    CHEP bad-channels=@{aux/files/bad-^.txt}  
  
  for (int p = 0 ; p < params.size(); p++ )
    params[p].update( id , globals::indiv_wildcard );
  
}


void cmd_t::replace_wildcards( const std::string & id )
{

  // ***** REDUNDANT *****
  // use the alternate form (defined above) instead of this legacy version
  //  - i.e. the code in this function can be removed when dynamic_parse()
  //         is well established
  // ***** REDUNDANT *****

  // ---->  i.e. use dynamic_parse() instead
  
  return dynamic_parse( id );


  // all text below if commented out now
  
  //
  // get copy of original script;  comments, newlines vesrus '&' will already have
  // been taken care of
  //
  
//  std::string iline = line;

  
  //
  // If the script contains an ID wildcard, confirm that ID does not already contain the wildcard 
  //
  
//  if ( iline.find( globals::indiv_wildcard ) != std::string::npos && 
//       id.find( globals::indiv_wildcard ) != std::string::npos ) 
//    Helper::halt( "ID " + id + " contains ID-wildcard character " 
//		  + globals::indiv_wildcard + " (i.e. use wildcard=X to specify in a different one)" );
 

  //
  // Copy a set of variables, where any i-vars will overwrite an existing var
  //
  
//  std::map<std::string,std::string> allvars = vars;
//  if ( ivars.find( id ) != ivars.end() )
//    {
//      const std::map<std::string,std::string> & newvars = ivars.find( id )->second;
//      std::map<std::string,std::string>::const_iterator vv = newvars.begin();      
//      while ( vv != newvars.end() )
//	{
//	  allvars[ vv->first ] = vv->second;
//	  ++vv;
//	}
//    }

  
  //
  // REMOVED:  IF/NOT/FI (or ENDIF) are now the canonical options 
  // remove any conditional blocks, where variable should be var=1 (in) or var=0 (out)
  // or a value for add=variable,var2 for var=1 var=2
  //
  
  // [[var1
  //   block
  // ]]var1  
  // Helper::process_block_conditionals( &iline , allvars );
  
  //
  // Now go line-by-line, to allow variables defined on-the-fly to feature in 
  // expansions, i.e. to build up expansions
  //

  // std::vector<std::string> cline = Helper::parse( iline , "\n" );
  
  // iline = "";
  
  // for (int l=0; l<cline.size(); l++)
  //   {

  //     std::string currline = cline[l];
      
  //     //
  //     // swap in any variables (and allows for them being defined on-the-fly)
  //     //
      
  //     Helper::swap_in_variables( &currline , &allvars );

  //     //
  //     // expand [term][1:10] [ == (term)[list] ] or [list](term) sequences
  //     //

  //     Helper::expand_numerics( &currline );
      
  //     iline += currline + "\n";
      
  //   }

  //
  // Parse into commands/options
  //

  //  std::vector<std::string> tok = Helper::quoted_parse( iline , "\n" );

  //
  // expand loops
  //

  //  tok = cmd_t::proc_forloops( tok );

  //
  // command(s): do this just for printing; real parsing will include variables, etc
  //
  
  // params.clear();
  // cmds.clear();

  // for (int c=0;c<tok.size();c++)
  //   {      
  //     std::vector<std::string> ctok = Helper::quoted_parse( tok[c] , "\t " );
      
  //     // may be 0 if a line was a variable declaration
  //     if ( ctok.size() >= 1 ) 
  // 	{
  // 	  //std::cout << "adding [" << ctok[0] << "]\n";
  // 	  cmds.push_back( ctok[0] );
	  
  // 	  param_t param;
  // 	  for (int j=1;j<ctok.size();j++) param.parse( ctok[j] );
  // 	  params.push_back( param );
  // 	}
  //   }
  
  // // replace in all 'params' any instances of 'globals::indiv_wildcard' with 'id'
  // // ALSO, will expand any @{includes} from files (which may contain ^ ID wildcards
  // // in their names,  e..g.    CHEP bad-channels=@{aux/files/bad-^.txt}  
  
  // for (int p = 0 ; p < params.size(); p++ )
  //   params[p].update( id , globals::indiv_wildcard );
   
}


bool cmd_t::read( const std::string * str , bool silent )
{
  
  bool cmdline_mode = str == NULL;   
  
  if ( std::cin.eof() && cmdline_mode ) return false;
  
  if ( (!cmdline_mode) && str->size() == 0 ) return false;
  
  reset();
  
  // CMD param=1 p=1,2,3 f=true out=o1 & CMD2 etc ; 
  
  // EITHER read from std::cin, 
  // OR          from -s command line
  // OR          from a string (e.g. lunaR)
  
  // Commands are delimited by & symbols (i.e. multi-line statements allowed)

  //   line continuations must start with + or ... (after whitespace)  
  // also, can have \ which means a continuation on the next line
  //  i.e replace \\n with space
  
  std::istringstream allinput;
  
  if ( ! cmdline_mode ) // read from 'str', such as R interface
    {
      
      // split by commands ('&' symbols), but allowing for & witin eval expressions 
      std::vector<std::string> tok = Helper::quoted_parse( *str , "&" );
      
      std::stringstream ss;
      
      for (int l=0;l<tok.size();l++)
	{
	  if ( tok[l] == "" ) continue;
	  if ( l != 0 ) ss << " & ";
	  ss << tok[l];
	}     
      allinput.str( ss.str() );

    }

  // read from std::cin
  else if ( cmd_t::cmdline_cmds == "" )
    {
      std::stringstream ss;
      bool first_cmd = true;
      while ( 1 )
	{
	  std::string s;
	  Helper::safe_getline( std::cin , s );
	  if ( std::cin.eof() ) break;
	  if ( s == "" ) continue;	  

	  s = Helper::ltrim( s ) ;
	  
	  // is this a continuation line?
	  // ... or + as first character
	  bool continuation = s.substr(0,1) == "+" || s.substr(0,3) == "..." ;
	  if ( continuation )
	    {
	      if ( s.substr(0,1) == "+" ) s = s.substr(1);
	      else s = s.substr(3);
	    }
	      
	  // only read up to a % comment, although this may be quoted
	  if ( s.find( "%" ) != std::string::npos ) 
	    {
	      bool inquote = false;
	      int comment_start = -1;
	      for (int i=0;i<s.size();i++)
		{
		  if ( s[i] == '"' ) inquote = ! inquote;
		  if ( s[i] == '%' && ! inquote ) { comment_start = i; break; }
		}
	      
	      // remove comment
	      if ( comment_start != -1 )
		s = s.substr( 0 , comment_start );
	    }

	  // trim leading/trailing whitespace
	  s = Helper::ltrim( s );
	  s = Helper::rtrim( s );

	  // anything left to add?
	  if ( s.size() > 0 ) 
	    {
	      if ( ! continuation ) 
		{
		  if ( ! first_cmd ) ss << " & ";		  
		  first_cmd = false;
		}
	      else 
		{
		  ss << " "; // spacer
		}	      
	      // add actual non-empty command
	      ss << s ;	      
	    }

	}

      allinput.str( ss.str() );
    }

  // read from -s string
  else
    {
      // allow multi-line commands via \ at end
      allinput.str( Helper::search_replace( cmd_t::cmdline_cmds , "\\\n", " " ) );
    }

  
  // take everything
  line = allinput.str();

  // change any '&' (back) to '\n', unless they are quoted (" or #)
  bool inquote = false;
  for (int i=0;i<line.size();i++) 
    {
      if ( line[i] == '"' ) inquote = ! inquote;
      else if ( line[i] == '&' ) 
	{
	  if ( ! inquote ) line[i] = '\n';
	}
    }

  // skip comments between commands if line starts with '%' (or '\n')
  bool recheck = true;
  while (recheck)
    {     
      if ( line[0] == '%' || line[0] == '\n' ) 
	{
	  line = line.substr( line.find("\n")+1);
	}
      else recheck = false;
    }
  

  //
  // Completed extracting 'line'
  //

  //
  // Note... we used to do processing of command file here (global),
  // but as we now use ivars as well as vars, hold on this till
  // later, in update()
  //
    

  //
  // Initial reporting of commands read
  //
  
  std::vector<std::string> tok = Helper::quoted_parse( line , "\n" );
  
  if ( tok.size() == 0 ) 
    {
      quit(true);      
      return false;
    }

    
  //
  // command(s): do this just for printing; real parsing will include variables, etc
  //
  
  for (int c=0;c<tok.size();c++)
    {      
      std::vector<std::string> ctok = Helper::quoted_parse( tok[c] , "\t " );
            
      // may be 0 if a line was a variable declaration
      if ( ctok.size() >= 1 ) 
	{
	  cmds.push_back( ctok[0] );
	  param_t param;
	  for (int j=1;j<ctok.size();j++) param.parse( ctok[j] );
	  params.push_back( param );
	}
    }

  
  // summary

  if ( ! silent )
    {
      logger << "input(s): " << input << "\n";
      logger << "output  : " << writer.name() 
	     << ( cmd_t::plaintext_mode ? " [dir for text-tables]" : "" ) 
	     << "\n";
      
      if ( signallist.size() > 0 )
	{
	  logger << "signals :";
	  for (std::set<std::string>::iterator s=signallist.begin();
	       s!=signallist.end();s++) 
	    logger << " " << *s ;
	  logger << "\n";
	}
      
      for (int i=0;i<cmds.size();i++)
	{
	  if ( i==0 ) 
	    logger << "commands: ";
	  else
	    logger << "        : ";
	  
	  logger << "c" << i+1 
		 << "\t" << cmds[i] << "\t"
		 << params[i].dump("","|") 
		 << "\n";
	}
    }

  return true;
}



//
// Evaluate commands
//

bool cmd_t::eval( edf_t & edf ) 
{
  
  //
  // Conditional evaluation
  //
  
  int if_count = 0;

  std::string if_condition = "";
  
  //
  // Loop over each command
  //
  
  for ( int c = 0 ; c < num_cmds() ; c++ )
    {	        
      
      // was a problem flag raised when loading the EDF?
          
      if ( globals::problem ) return false;
      
      //
      // If this particular command did not explicitly specify
      // signals, then add a wildcard, which will always take all
      // in-memory channels
      //
            
      if ( ! param(c).has( "sig" ) )
	param(c).add_hidden( "sig" , "*" );
      
      //
      // Print command
      //

      logger << " ..................................................................\n"
	     << " CMD #" << c+1 << ": " << cmd(c) << "\n";
      logger << "   options: " << param(c).dump( "" , " " ) << "\n";

      
      //
      // Quit?
      //

      if ( is( c, "EXIT" ) || 
	   is( c, "STOP" ) || 
	   is( c, "QUIT" ) ) 
	{
	  param_t p = param(c);
	  const bool req_problem = p.has( "if-problem" ) || p.has( "problem" );
	  const bool req_empty = p.has( "if-empty" ) || p.has( "empty" );
	  	  	  
	  // just quit?
	  if ( ! ( req_problem || req_empty ) )
	    return true;

	  // else, only quit if we match a problem
	  if ( req_problem && globals::problem ) return true;
	  if ( req_empty && globals::empty ) return true;
	  
	  // otherwise, carry on
	  logger << "  not quitting as no problem flag has been set\n";
	  continue;
	}


      //
      // Is the current mask empty? if so, skip unless this is a THAW
      //

      if ( globals::empty )
	{
	  // not systematic, but some commands that do not access the EDF are okay 
	  // to rerun - most importantly THAW... 
	  bool skip = true;
	  if      ( is( c, "THAW" ) ) skip = false;
	  else if ( is( c, "TAG" ) ) skip = false; 
	  else if ( is( c, "HEADERS" ) ) skip = false;
	  else if ( is( c, "SET-VAR" ) ) skip = false;
	  else if ( is( c, "SET-HEADERS" ) ) skip = false;
	  else if ( is( c, "DESC" ) ) skip = false;
	  else if ( is( c, "ALIASES" ) ) skip = false;
	  else if ( is( c, "TYPES" ) ) skip = false;
	  else if ( is( c, "PREDICT" ) ) skip = false;
	  else if ( is( c, "REPORT" ) ) skip = false;
	  if ( skip )
	    {
	      logger << "  ** skipping " << cmd(c) << " as there are no unmasked records\n"; 
	      continue;
	    }
	}


      //
      // Check that epoch type is consistent w/ the command
      //
      
      if ( ! edf.timeline.check( cmd(c) ) )
	Helper::halt( "cannot apply " + cmd(c) + " with non-standard epoch definitions" );
      
      //
      // Process command
      //
      
      writer.cmd( cmd(c) , c+1 , param(c).dump( "" , " " ) );
      
      // use strata to keep track of tables by commands, with leading underscore to denote
      // special status of this factor
      
      writer.level( cmd(c) , "_" + cmd(c) );

      
      //
      // Now process the command
      //
      bool fnd = false;
      
      if ( (!fnd) && is( c, "WRITE" ) )        { fnd = true; proc_write( edf, param(c) ); }
      if ( (!fnd) && is( c, "EDF-" ) )         { fnd = true; proc_edf_minus( edf , param(c) ); }
      if ( (!fnd) && is( c, "EDF-MINUS" ) )    { fnd = true; proc_edf_minus( edf , param(c) ); }
      if ( (!fnd) && is( c, "SET-TIMESTAMPS")) { fnd = true; proc_set_timestamps( edf , param(c) ); }
      if ( (!fnd) && is( c, "DROP-TIMESTAMPS")){ fnd = true; proc_force_edf( edf , param(c) ); }      
      if ( (!fnd) && is( c, "SUMMARY" ) )      { fnd = true; proc_summaries( edf , param(c) ); }
      if ( (!fnd) && is( c, "HEADERS" ) )      { fnd = true; proc_headers( edf , param(c) ); }
      if ( (!fnd) && is( c, "ALIASES" ) )      { fnd = true; proc_aliases( edf , param(c) ); }
      if ( (!fnd) && is( c, "REPORT" ) )       { fnd = true; proc_report( edf , param(c) ); }
      if ( (!fnd) && is( c, "SET-HEADERS" ) )  { fnd = true; proc_set_headers( edf , param(c) ); }
      if ( (!fnd) && is( c, "SET-VAR" ) )      { fnd = true; proc_set_ivar( edf, param(c) ); }
      if ( (!fnd) && is( c, "DESC" ) )         { fnd = true; proc_desc( edf , param(c) ); }
      if ( (!fnd) && is( c, "TYPES" ) )        { fnd = true; proc_show_channel_map(); }
      if ( (!fnd) && is( c, "VARS" ) )         { fnd = true; proc_dump_vars( edf , param(c) ); }
      if ( (!fnd) && is( c, "STATS" ) )        { fnd = true; proc_stats( edf , param(c) ); }
      if ( (!fnd) && is( c, "DUPES" ) )        { fnd = true; proc_dupes( edf, param(c) );  }
      if ( (!fnd) && is( c, "REFERENCE" ) )    { fnd = true; proc_reference( edf , param(c) ); }
      if ( (!fnd) && is( c, "DEREFERENCE" ) )  { fnd = true; proc_dereference( edf , param(c) ); }      
      if ( (!fnd) && is( c, "FLIP" ) )         { fnd = true; proc_flip( edf , param(c) ); }
      if ( (!fnd) && is( c, "RAI" ) )          { fnd = true; proc_rai( edf, param(c) ); }
      if ( (!fnd) && is( c, "RECTIFY" ) )      { fnd = true; proc_rectify( edf , param(c) ); }
      if ( (!fnd) && is( c, "CLIP" ) )         { fnd = true; proc_clip( edf, param(c) ); } 
      if ( (!fnd) && is( c, "REVERSE" ) )      { fnd = true; proc_reverse( edf , param(c) ); }
      if ( (!fnd) && is( c, "CANONICAL" ) )    { fnd = true; proc_canonical( edf , param(c) ); }
      if ( (!fnd) && is( c, "REMAP" ) )        { fnd = true; proc_remap_annots( edf , param(c) ); }
      if ( (!fnd) && is( c, "uV" ) )           { fnd = true; proc_scale( edf , param(c) , "uV" ); }
      if ( (!fnd) && is( c, "mV" ) )           { fnd = true; proc_scale( edf , param(c) , "mV" ); }
      if ( (!fnd) && is( c, "MINMAX" ) )       { fnd = true; proc_minmax( edf , param(c) ); }
      if ( (!fnd) && is( c, "SCALE" ) )        { fnd = true; proc_setscale( edf , param(c) ); }
      if ( (!fnd) && is( c, "ROBUST-NORM" ) )  { fnd = true; proc_standardize( edf , param(c) ); }
      if ( (!fnd) && is( c, "ROLLING-NORM" ) ) { fnd = true; proc_rolling_norm( edf, param(c) ); }
      if ( (!fnd) && is( c, "ALTER" ) )        { fnd = true; proc_correct( edf , param(c) ); }
      if ( (!fnd) && is( c, "RECORD-SIZE" ) )  { fnd = true; proc_rerecord( edf , param(c) ); }  
      if ( (!fnd) && is( c, "TIME-TRACK" ) )   { fnd = true; proc_timetrack( edf, param(c) ); }
      if ( (!fnd) && is( c, "STAGE" ) )        { fnd = true; proc_sleep_stage( edf , param(c) , false ); }
      if ( (!fnd) && is( c, "HYPNO" ) )        { fnd = true; proc_sleep_stage( edf , param(c) , true ); }
      if ( (!fnd) && is( c, "SIG-CYCLES") )    { fnd = true; proc_ecycle( edf , param(c) ); }
      if ( (!fnd) && is( c, "TSLIB" ) )        { fnd = true; pdc_t::construct_tslib( edf , param(c) ); }
      if ( (!fnd) && is( c, "SSS" ) )          { fnd = true; pdc_t::simple_sleep_scorer( edf , param(c) ); }
      if ( (!fnd) && is( c, "EXE" ) )          { fnd = true; pdc_t::similarity_matrix( edf , param(c) ); }
      if ( (!fnd) && is( c, "DUMP" ) )         { fnd = true; proc_dump( edf, param(c) );  }
      if ( (!fnd) && is( c, "DUMP-RECORDS" ) ) { fnd = true; proc_record_dump( edf , param(c) ); }
      if ( (!fnd) && is( c, "RECS" ) )         { fnd = true; proc_record_table( edf , param(c) ); }
      if ( (!fnd) && is( c, "SEGMENTS" ) )     { fnd = true; proc_dump_segs( edf , param(c) ); }
      if ( (!fnd) && is( c, "DUMP-EPOCHS" ) )  { fnd = true; proc_epoch_dump( edf, param(c) ); } // REDUNDANT -> ANNOTS epoch
      if ( (!fnd) && is( c, "ANNOTS" ) )       { fnd = true; proc_list_all_annots( edf, param(c) ); }
      if ( (!fnd) && is( c, "ESPAN" ) )        { fnd = true; proc_espan( edf , param(c) ); }
      if ( (!fnd) && is( c, "MAKE-ANNOTS" ) )  { fnd = true; proc_make_annots( edf , param(c) ); }
      if ( (!fnd) && is( c, "WRITE-ANNOTS" ) ) { fnd = true; proc_write_annots( edf, param(c) ); }
      if ( (!fnd) && is( c, "META" ) )         { fnd = true; proc_set_annot_metadata( edf, param(c) );  }
      if ( (!fnd) && is( c, "OVERLAP") )       { fnd = true; proc_annotate( edf, param(c) ); }
      if ( (!fnd) && is( c, "EXTEND" ) )       { fnd = true; proc_extend_annots( edf, param(c) ); }
      if ( (!fnd) && is( c, "AXA") )           { fnd = true; proc_annot_crosstabs( edf, param(c) ); }
      if ( (!fnd) && is( c, "A2S" ) )          { fnd = true; proc_annot2signal( edf, param(c) ); }
      if ( (!fnd) && is( c, "S2A" ) )          { fnd = true; proc_signal2annot( edf, param(c) ); }
      if ( (!fnd) && is( c, "A2C" ) )          { fnd = true; proc_annot2cache( edf , param(c) ); }
      if ( (!fnd) && is( c, "C2A" ) )          { fnd = true; proc_cache2annot( edf , param(c) ); }
      if ( (!fnd) && is( c, "SPANNING" ) )     { fnd = true; proc_list_spanning_annots( edf, param(c) ); }
      if ( (!fnd) && is( c, "MEANS" ) )        { fnd = true; proc_sig_annot_mean( edf, param(c) ); }
      if ( (!fnd) && is( c, "TABULATE" ) )     { fnd = true; proc_sig_tabulate( edf, param(c) ); }
      if ( (!fnd) && is( c, "MATRIX" ) )       { fnd = true; proc_epoch_matrix( edf , param(c) ); }
      if ( (!fnd) && is( c, "HEAD" ) )         { fnd = true; proc_head_matrix( edf , param(c) ); }
      if ( (!fnd) && ( is( c, "RESTRUCTURE" ) || is( c, "RE" ) ) )  { fnd = true; proc_restructure( edf , param(c) ); }
      if ( (!fnd) && is( c, "SIGNALS" ) )      { fnd = true; proc_drop_signals( edf , param(c) ); }
      if ( (!fnd) && is( c, "DROP-ANNOTS" ) )  { fnd = true; proc_drop_annots( edf , param(c) ); }
      if ( (!fnd) && is( c, "RENAME" ) )       { fnd = true; proc_rename( edf, param(c) ); }
      if ( (!fnd) && is( c, "ENFORCE-SR" ) )   { fnd = true; proc_enforce_signals( edf , param(c) ); }
      if ( (!fnd) && is( c, "COPY" ) )         { fnd = true; proc_copy_signal( edf , param(c) ); }
      if ( (!fnd) && is( c, "READ" ) )         { fnd = true; proc_read_signal( edf , param(c) ); }
      if ( (!fnd) && is( c, "ORDER" ) )        { fnd = true; proc_order_signals( edf , param(c) ); }
      if ( (!fnd) && is( c, "REQUIRES" ) )     { fnd = true; proc_requires( edf , param(c) ); } 
      if ( (!fnd) && is( c, "CONTAINS" ) )     { fnd = true; proc_has_signals( edf , param(c) ); }
      if ( (!fnd) && ( is( c, "RMS" ) || is( c, "SIGSTATS" ) ) ) { fnd = true; proc_rms( edf, param(c) ); }
      if ( (!fnd) && is( c, "MSE" ) )          { fnd = true; proc_mse( edf, param(c) ); }
      if ( (!fnd) && is( c, "LZW" ) )          { fnd = true; proc_lzw( edf, param(c) ); }
      if ( (!fnd) && is( c, "ZR" ) )           { fnd = true; proc_zratio( edf , param(c) ); }
      if ( (!fnd) && is( c, "ANON" ) )         { fnd = true; proc_anon( edf , param(c) ); }
      if ( (!fnd) && is( c, "EPOCH" ) )        { fnd = true; proc_epoch( edf, param(c) ); }
      if ( (!fnd) && is( c, "ALIGN" ) )        { fnd = true; proc_align( edf , param(c) ); }
      if ( (!fnd) && is( c, "SLICE" ) )        { fnd = true; proc_slice( edf , param(c) , 1 ); }
      if ( (!fnd) && is( c, "ALIGN-EPOCHS" ) ) { fnd = true; proc_align_epochs( edf , param(c) ); }
      if ( (!fnd) && is( c, "ALIGN-ANNOTS" ) ) { fnd = true; proc_align_annots( edf , param(c) ); }
      if ( (!fnd) && is( c, "INSERT" ) )       { fnd = true; proc_insert( edf , param(c) ); }
      if ( (!fnd) && is( c, "SUDS" ) )         { fnd = true; proc_suds( edf , param(c) ); }
      if ( (!fnd) && is( c, "MAKE-SUDS" ) )    { fnd = true; proc_make_suds( edf , param(c) ); }
      if ( (!fnd) && is( c, "RUN-POPS" ) )     { fnd = true; proc_runpops( edf , param(c) ); } // wrapper
      if ( (!fnd) && is( c, "POPS" ) )         { fnd = true; proc_pops( edf , param(c) ); }
      if ( (!fnd) && is( c, "EVAL-STAGES" ) )  { fnd = true; proc_eval_stages( edf , param(c) ); }
      if ( (!fnd) && is( c, "SOAP" ) )         { fnd = true; proc_self_suds( edf , param(c) ); }
      if ( (!fnd) && is( c, "COMPLETE" ) )     { fnd = true; proc_resoap( edf , param(c) ); }
      if ( (!fnd) && is( c, "REBASE" ) )       { fnd = true; proc_rebase_soap( edf , param(c) ); } // e.g. 20->30s epochs using SOAP
      if ( (!fnd) && is( c, "PLACE" ) )        { fnd = true; proc_place_soap( edf , param(c) ); } // e.g. find where should go
      if ( (!fnd) && is( c, "TRANS" ) )        { fnd = true; proc_trans( edf , param(c) ); }
      if ( (!fnd) && is( c, "EVAL" ) )         { fnd = true; proc_eval( edf, param(c) ); }
      if ( (!fnd) && is( c, "DERIVE" ) )       { fnd = true; proc_derive( edf, param(c) ); }
      if ( (!fnd) && is( c, "MASK" ) )         { fnd = true; proc_mask( edf, param(c) ); }
      if ( (!fnd) && is( c, "COMBINE" ) )      { fnd = true; proc_combine( edf , param(c) ); }
      if ( (!fnd) && is( c, "FREEZE" ) )       { fnd = true; proc_freeze( edf , param(c) ); }
      if ( (!fnd) && is( c, "THAW" ) )         { fnd = true; proc_thaw( edf , param(c) ); }
      if ( (!fnd) && is( c, "CLEAN-FREEZER"))  { fnd = true; proc_clean_freezer( edf , param(c) ); }
      if ( (!fnd) && is( c, "FILE-MASK" ) )    { fnd = true; proc_file_mask( edf , param(c) ); } // not supported/implemented
      if ( (!fnd) && is( c, "DUMP-MASK" ) )    { fnd = true; proc_dump_mask( edf, param(c) ); }
      if ( (!fnd) && is( c, "ANNOT-MASK" ) )   { fnd = true; proc_annot_mask( edf, param(c) ); }
      if ( (!fnd) && is( c, "CHEP" ) )         { fnd = true; timeline_t::proc_chep( edf, param(c) ); }
      if ( (!fnd) && is( c, "CHEP-MASK" ) )    { fnd = true; proc_chep_mask( edf, param(c) ); }
      if ( (!fnd) && is( c, "EPOCH-ANNOT" ) )  { fnd = true; proc_file_annot( edf , param(c) ); }
      if ( (!fnd) && is( c, "EPOCH-MASK" ) )   { fnd = true; proc_epoch_mask( edf, param(c) ); }
      if ( (!fnd) && is( c, "HB" ) )           { fnd = true; proc_hypoxic_burden( edf, param(c) ); }
      if ( (!fnd) && is( c, "FILTER" ) )       { fnd = true; proc_filter( edf, param(c) ); }
      if ( (!fnd) && is( c, "FILTER-DESIGN" )) { fnd = true; proc_filter_design( edf, param(c) ); }
      if ( (!fnd) && is( c, "MOVING-AVERAGE" )) { fnd = true; proc_moving_average( edf, param(c) ); }
      if ( (!fnd) && is( c, "CWT-DESIGN" ) )   { fnd = true; proc_cwt_design( edf , param(c) ); }
      if ( (!fnd) && is( c, "CWT" ) )          { fnd = true; proc_cwt( edf , param(c) ); }
      if ( (!fnd) && is( c, "HILBERT" ) )      { fnd = true; proc_hilbert( edf , param(c) ); }
      if ( (!fnd) && is( c, "SYNC" ) )         { fnd = true; proc_sync( edf , param(c) ); }
      if ( (!fnd) && is( c, "TSYNC" ) )        { fnd = true; proc_tsync( edf , param(c) ); }
      if ( (!fnd) && is( c, "XCORR" ) )        { fnd = true; proc_xcorr( edf, param(c) ); }
      if ( (!fnd) && is( c, "TV" ) )           { fnd = true; proc_tv_denoise( edf , param(c) ); }
      if ( (!fnd) && is( c, "OTSU" ) )         { fnd = true; proc_otsu( edf, param(c) ); }
      if ( (!fnd) && is( c, "COVAR" ) )        { fnd = true; proc_covar( edf, param(c) ); }
      if ( (!fnd) && is( c, "PSD" ) )          { fnd = true; proc_psd( edf, param(c) ); }
      if ( (!fnd) && is( c, "FFT" ) )          { fnd = true; proc_fft( edf , param(c) ); }
      if ( (!fnd) && is( c, "MTM" ) )          { fnd = true; proc_mtm( edf, param(c) ); }
      if ( (!fnd) && is( c, "IRASA" ) )        { fnd = true; proc_irasa( edf, param(c) ); }
      if ( (!fnd) && is( c, "1FNORM" ) )       { fnd = true; proc_1overf_norm( edf, param(c) ); }
      if ( (!fnd) && is( c, "DYNAM" ) )        { fnd = true; proc_qdynam( edf , param(c) ); }
      if ( (!fnd) && is( c, "PSC" ) )          { fnd = true; proc_psc( edf , param(c) ); }
      if ( (!fnd) && is( c, "MS" ) )           { fnd = true; proc_microstates( edf , param(c) ); }
      if ( (!fnd) && is( c, "ASYMM" ) )        { fnd = true; proc_asymm( edf , param(c) ); }
      if ( (!fnd) && is( c, "TLOCK" ) )        { fnd = true; proc_tlock( edf , param(c) ); }
      if ( (!fnd) && is( c, "TCLST" ) )        { fnd = true; proc_tclst( edf , param(c) ); }
      if ( (!fnd) && is( c, "PERI" ) )         { fnd = true; proc_peri( edf , param(c) ); }
      if ( (!fnd) && is( c, "PEAKS" ) )        { fnd = true; proc_peaks( edf , param(c) ); } 
      if ( (!fnd) && is( c, "Z-PEAKS" ) )      { fnd = true; proc_zpeaks( edf , param(c) ); }
      if ( (!fnd) && is( c, "SEDF" ) )         { fnd = true; proc_sedf( edf , param(c) ); }
      if ( (!fnd) && is( c, "FIP" ) )          { fnd = true; proc_fiplot( edf , param(c) ); }
      if ( (!fnd) && is( c, "COH" ) )          { fnd = true; proc_coh( edf , param(c) ); }
      if ( (!fnd) && is( c, "CC" ) )           { fnd = true; proc_conncoupl( edf , param(c) ); }
      if ( (!fnd) && is( c, "CORREL" ) )       { fnd = true; proc_correl( edf , param(c) ); }
      if ( (!fnd) && is( c, "PSI" ) )          { fnd = true; proc_psi( edf , param(c) ); }
      if ( (!fnd) && is( c, "ACF" ) )          { fnd = true; proc_acf( edf , param(c) ); }
      if ( (!fnd) && is( c, "GP" ) )           { fnd = true; gc_wrapper( edf , param(c) );  }
      if ( (!fnd) && is( c, "ED" ) )           { fnd = true; proc_elec_distance( edf , param(c) ); }
      if ( (!fnd) && is( c, "SVD" ) )          { fnd = true; proc_svd( edf, param(c) ); }
      if ( (!fnd) && is( c, "ICA" ) )          { fnd = true; proc_ica( edf, param(c) ); }
      if ( (!fnd) && is( c, "ADJUST" ) )       { fnd = true; proc_adjust( edf , param(c) );  }
      if ( (!fnd) && is( c, "CLOCS" ) )        { fnd = true; proc_attach_clocs( edf , param(c) ); }
      if ( (!fnd) && is( c, "L1OUT" ) )        { fnd = true; proc_leave_one_out( edf , param(c) ); }
      if ( (!fnd) && is( c, "INTERPOLATE" ) )  { fnd = true; proc_chep_based_interpolation( edf, param(c) ); }
      if ( (!fnd) && is( c, "SL" ) )           { fnd = true; proc_surface_laplacian( edf , param(c) ); }
      if ( (!fnd) && is( c, "EMD" ) )          { fnd = true; proc_emd( edf , param(c) ); }
      if ( (!fnd) && is( c, "DFA" ) )          { fnd = true; proc_dfa( edf , param(c) ); }
      if ( (!fnd) && is( c, "MI" ) )           { fnd = true; proc_mi( edf, param(c) ); }
      if ( (!fnd) && is( c, "HR" ) )           { fnd = true; proc_bpm( edf , param(c) ); }
      if ( (!fnd) && is( c, "HRV" ) )          { fnd = true; proc_hrv( edf , param(c) ); }
      if ( (!fnd) && is( c, "SUPPRESS-ECG" ) ) { fnd = true; proc_ecgsuppression( edf , param(c) ); }
      if ( (!fnd) && is( c, "PAC" ) )          { fnd = true; proc_pac( edf , param(c) ); }
      if ( (!fnd) && is( c, "CFC" ) )          { fnd = true; proc_cfc( edf , param(c) ); }
      if ( (!fnd) && is( c, "GED" ) )          { fnd = true; proc_ged( edf , param(c) ); }
      if ( (!fnd) && is( c, "PREDICT" ) )      { fnd = true; proc_predict( edf , param(c) ); }
      if ( (!fnd) && is( c, "TAG" ) )          { fnd = true; proc_tag( param(c) ); }
      if ( (!fnd) && is( c, "RESAMPLE" ) )     { fnd = true; proc_resample( edf, param(c) ); }
      if ( (!fnd) && is( c, "FIX-SAMPLING" ) ) { fnd = true; proc_fix_sampling( edf , param(c) ); }
      if ( (!fnd) && is( c, "ZOH" ) )          { fnd = true; proc_zoh( edf, param(c) ); }
      if ( (!fnd) && is( c, "LINE-DENOISE" ) ) { fnd = true; dsptools::line_denoiser( edf, param(c) ); }
      if ( (!fnd) && is( c, "ZC" ) )           { fnd = true; dsptools::detrend( edf, param(c) ); }
      if ( (!fnd) && is( c, "SPINDLES" ) )     { fnd = true; proc_spindles( edf, param(c) ); }
      if ( (!fnd) && is( c, "AROUSALS" ) )     { fnd = true; proc_arousals( edf, param(c) ); } 
      if ( (!fnd) && is( c, "SO" ) )           { fnd = true; proc_slowwaves( edf, param(c) ); }
      if ( (!fnd) && is( c, "COUPL" ) )        { fnd = true; proc_coupling( edf , param(c) ); }
      if ( (!fnd) && is( c, "RIPPLES" ) )      { fnd = true; proc_ripples( edf , param(c) ); }
      if ( (!fnd) && is( c, "PCOUPL" ) )       { fnd = true; proc_generic_coupling( edf , param(c) ); }
      if ( (!fnd) && is( c, "POL" ) )          { fnd = true; proc_polarity( edf, param(c) ); }	  
      if ( (!fnd) && is( c, "REMS" ) )         { fnd = true; proc_rems( edf, param(c) ); }
      if ( (!fnd) && is( c, "QC" ) )           { fnd = true; proc_qc( edf, param(c) ); }
      if ( (!fnd) && is( c, "ARTIFACTS" ) )    { fnd = true; proc_artifacts( edf, param(c) ); }
      if ( (!fnd) && is( c, "EDGER" ) )        { fnd = true; proc_trim( edf, param(c) ) ; }
      if ( (!fnd) && is( c, "CACHE" ) )        { fnd = true; proc_dump_cache( edf , param(c) ); }
      if ( (!fnd) && is( c, "SIGGEN" ) )       { fnd = true; proc_siggen( edf , param(c) ); }
      if ( (!fnd) && is( c, "SIMUL" ) )        { fnd = true; proc_simul( edf , param(c) ); }
      if ( (!fnd) && is( c, "SPIKE" ) )        { fnd = true; proc_spike( edf , param(c) ); }
      if ( (!fnd) && is( c, "SHIFT" ) )        { fnd = true; proc_shift( edf , param(c) ); }
      if ( (!fnd) && is( c, "SCRAMBLE" ) )     { fnd = true; proc_scramble( edf , param(c) ); }
      
#ifdef HAS_LGBM
      if ( (!fnd) && is( c, "PREP-MASSOC" ) )  { fnd = true; massoc_t::massoc_dumper( edf , param(c) ); }
#endif
      
      //
      // we don't recognize this command
      //
      
      if ( ! fnd ) 
	{
	  Helper::halt( "did not recognize command: " + cmd(c) );
	  return false; 
	}

       
      //
      // Was a problem flag set?
      //
    
      if ( globals::problem ) 
	{
	  
	  logger << "**warning: the PROBLEM flag was set, skipping to next EDF...\n";

	  writer.value( "PROBLEM" , 1 );

	  if ( globals::write_naughty_list )
	    {
	      logger << "**writing ID " << edf.id << " to " <<  globals::naughty_list << "\n";
	      std::ofstream PROBLEMS( globals::naughty_list.c_str() , std::ios_base::app );
	      PROBLEMS << edf.id << "\n";
	      PROBLEMS.close();
	    }
	  
	  return false;
	}
      
      //
      // next command
      //

      writer.unlevel( "_" + cmd(c) );      
     

    } // next command
  


  return true;
}



//
// ----------------------------------------------------------------------------------------
//
// Wrapper (proc_*) functions below that drive specific commands
//
//  - most simply pass edf and param to the handling function
//  - in some (generally older) cases, there is some pre-processing done here
//  - these are not ordered in any type of logical manner...
//
// ----------------------------------------------------------------------------------------
//

// HEADERS : summarize EDF files
void proc_headers( edf_t & edf , param_t & param )
{
  // optionally add a SIGNALS col that has a comma-delimited list of all signals
  edf.terse_summary( param );
}

// SET-VAR : set an IVAR
void proc_set_ivar( edf_t & edf , param_t & param )
{
  // expect var=" expr "
  
  std::string expr;

  std::string var = param.single_pair( & expr );
  
  logger << "  setting individual-level variable " << var << " to [ " << expr << " ]\n";
  
  // read a single line
  std::map<std::string,annot_map_t> inputs;
  
  instance_t out;
  
  Eval tok( expr );

  tok.bind( inputs , &out );

  const bool verbose = true;
  
  bool is_valid = tok.evaluate( verbose );
  
  bool retval;
  
  if ( ! tok.value( retval ) ) is_valid = false;

  std::cout << "parsed as a valid expression : " << ( is_valid ? "yes" : "no" ) << "\n";
  std::cout << "return value                 : " << tok.result() << "\n";
  std::cout << "return value (as T/F)        : " << ( retval ? "true" : "false" ) << "\n";
  std::cout << "assigned meta-data           : " << out.print() << "\n";  

  // and set (as string value)
  if ( is_valid ) 
    cmd_t::ivars[ edf.id ][ var ] = tok.result();
}

// SET-HEADERS : set EDF header fields
void proc_set_headers( edf_t & edf , param_t & param )
{
  edf.set_headers( param );
}

// ALIASES : report aliasing of channels and annotations
void proc_aliases( edf_t & edf , param_t & param )
{
  edf.report_aliases();
}
    
// REPORT : ensure VAR is added in output
void proc_report( edf_t & edf , param_t & param )
{

  // two functions::
  //  1) ensure a variable/ftable is added to the output, so that text-table works (ensure)
  //  2) sets hide/show for a cmd/fac/var/all to hide/show from all output streams (hide/show)

  // the first (ensure) is really only a kludge to handle cases where the table/variable
  // has not been added to cmddefs.
  
  const bool hide_mode = param.has( "hide" );
  const bool show_mode = param.has( "show" );

  const bool hide_all = param.has( "hide-all" );

  const bool show_all = param.has( "show-all" );
  
  const bool ensure_mode = param.has( "ensure" );

  const bool compress_mode = param.has( "compress" );
  
  
  if ( ! ( compress_mode || ensure_mode || hide_mode || show_mode || hide_all || show_all ) )
    Helper::halt( "nothing to do with REPORT: specify hide/show, hide-all/show-all, compress or ensure" );

  if ( ( compress_mode || ensure_mode ) && ( hide_mode || show_mode || hide_all || show_all ) )
    Helper::halt( "if using ensure mode, cannot specify hide/show commands" );

  if ( compress_mode && ensure_mode )
    Helper::halt( "cannot specify both compress and ensure" );
  
  
  // allow hide-all or show-all to be combined w/ a show or hide command

  if ( hide_mode && show_mode )
    Helper::halt( "can only specify one of hide and show" );
  
  if ( hide_all && show_all )
    Helper::halt( "can only specify one of hide-all and show-all" );
  
  if ( hide_all )
    {
      logger << "  hiding all outputs\n";
      globals::cmddefs().hide_all();
    }

  if ( show_all )
    {
      logger << "  showing all outputs\n";
      globals::cmddefs().show_all();
    }
  
  // at this point, will always require cmd,
  // optionally, fac
  // optionally vars
  
  const std::string cmd = param.has( "cmd" ) ? param.value( "cmd" ) : "";
    
  const bool has_table = param.has( "fac" );

  const bool has_vars = param.has( "vars" );
  
  if ( ( ! ( hide_all || show_all ) ) || has_table || has_vars )
    if ( ! param.has( "cmd" ) )
      Helper::halt( "no 'cmd' argument specified" );
  
  if ( has_vars && ! has_table )
    Helper::halt( "vars specified but no fac" );
  
  std::string fac;
  if ( has_table )
    {
      fac = param.value( "fac" );
      // allow BL or "." to specify baseline ( i.e. "" ) in cmddefs.
      if ( fac == "BL" || fac == "." || fac == "bl" ) fac = "";
    }
  
  std::vector<std::string> vars;
  if ( has_vars )
    vars = param.strvector( "vars" );

  //
  // compress mode
  //

  if ( compress_mode )
    {
      if ( has_table )
	{
	  logger << "  setting " << cmd << "/" << fac << " text-table output to compressed\n";
	  globals::cmddefs().set_compressed( cmd , tfac_t( fac ) );
	}
      else
        {
	  logger << "  setting " << cmd << " text-table output to compressed\n";
	  globals::cmddefs().set_compressed( cmd );
	}
      return;
    }
  
  //
  // logging check
  //
  
  if ( ensure_mode )
    {
      for (int v=0; v<vars.size(); v++)
	{
	  logger << "  logging " << vars[v]
		 << " for " << cmd << ( fac != "" ? ": " + fac : "" ) << "\n";
	  
	  globals::cmddefs().ensure_table( cmd , fac );
	  globals::cmddefs().add_var( cmd , fac , vars[v] , "." );
	}
      return;
    }

  //
  // else, we are setting either cmd, a table or specify variable(s) to show or not
  //
  
  if ( hide_mode )
    {      
      if ( has_vars )
	for (int v=0; v<vars.size(); v++)
	  globals::cmddefs().hide_var( cmd , fac , vars[v] );	
      else if ( has_table )
	globals::cmddefs().hide_table( cmd , fac );	
      else
	globals::cmddefs().hide_cmd( cmd );

      return;
    }


  if ( show_mode )
    {      
      if ( has_vars )
	for (int v=0; v<vars.size(); v++)
	  globals::cmddefs().show_var( cmd , fac , vars[v] );	
      else if ( has_table )
	globals::cmddefs().show_table( cmd , fac );	
      else
	globals::cmddefs().show_cmd( cmd );

      return;
    }

}


// SUMMARY : summarize EDF files (verbose, human-readable)  
void proc_summaries( edf_t & edf , param_t & param )
{
  std::cout << "EDF filename   : " << edf.filename << "\n" 
	    << edf.header.summary() << "\n"
	    << "----------------------------------------------------------------\n\n";
}


// DESC : very brief summary of contents 
void proc_desc( edf_t & edf , param_t & param )
{
  edf.description( param );
}

// TYPES : show channel mappings
void proc_show_channel_map()
{
  std::cout << globals::dump_channel_map(); 
}

// VARS : dump variables for this individual
void proc_dump_vars( edf_t & edf , param_t & param )
{
  
  std::map<std::string,std::string>::const_iterator vv =  cmd_t::vars.begin();
  while ( vv != cmd_t::vars.end() )
    {
      writer.level( vv->first , globals::var_strat );
      writer.value( "INDIV" , 0 );
      writer.value( "VAL" , vv->second );
      ++vv;
    }
  writer.unlevel( globals::var_strat );

  // no individual-level variables
  if ( cmd_t::ivars.find( edf.id ) == cmd_t::ivars.end() ) return;
  
  const std::map<std::string,std::string> & ivars =  cmd_t::ivars.find( edf.id )->second;

  vv = ivars.begin();
  while ( vv != ivars.end() )
    {
      writer.level( vv->first , globals::var_strat );
      writer.value( "INDIV" , 1 );
      writer.value( "VAL" , vv->second );
      ++vv;
    }
  writer.unlevel( globals::var_strat );
        
}

// DUPES : find signals that are ~duplicates
void proc_dupes( edf_t & edf , param_t & param )
{
  dsptools::dupes( edf, param );
}

// STATS : get basic stats for an EDF
void proc_stats( edf_t & edf , param_t & param )
{
  edf.basic_stats( param );
}

// RMS/SIGSTATS : calculate root-mean square for signals 
void proc_rms( edf_t & edf , param_t & param )
{    
  rms_per_epoch( edf , param );
}

// MSE : calculate multi-scale entropy per epoch
void proc_mse( edf_t & edf , param_t & param )
{    
  mse_per_epoch( edf , param );
}

// LZW : compression index per epoch/signal
void proc_lzw( edf_t & edf , param_t & param )
{    
  lzw_per_epoch( edf , param );
}


// SOAP : single observation accuracies and probabilities;
//   i.e. use SUDS on self, to evaluate staging/signal quality
void proc_self_suds( edf_t & edf , param_t & param  )
{

  // two modes:
  //    - if using the default model, and multiple signals, run once
  //      for each
  //    - otherwise, force signal 

  // if multiple channels specified, run once for each
  signal_list_t signals = edf.header.signal_list( param.value( "sig" ) );
  
  if ( signals.size() == 0 )
    Helper::halt( "no signals found matching " + param.value( "sig" ) );

  // force a reload (used in moonlight/R mode)
  if ( param.has( "force-reload" ) )
    suds_t::model.init();
  
  // set options
  suds_t::set_options( param );  

  const int ns = signals.size();

  for (int s=0; s<ns; s++)
    {
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;

      logger << "  ------------------------------------------------------------\n"
	     << "  fitting SOAP model for single channel " << signals.label(s) << "\n";
      
      writer.level( signals.label(s) , globals::signal_strat );
      
      // force reload either way
      suds_t::model.init();
            
      suds_t::model.read( param.has( "model" ) ? param.value( "model" ) : "_1" , 
			  param.has( "read-weights" ) ? param.value( "read-weights" ) : "" ,
			  param.has( "write-weights" ) ? param.value( "write-weights" ) : "" ,  
			  signals.label(s) );
      
      suds_indiv_t self;

      self.evaluate( edf , param );
      
      // next channel
    }

  writer.unlevel( globals::signal_strat );
    
}


// PLACE : find where stages should go
void proc_place_soap( edf_t & edf , param_t & param  )
{
  // expected .eannot style file
  std::string stagefile = param.requires( "stages" );  

  suds_t::set_options( param );

  // load model, if not already done (or default) (_1 or _2)

  if ( ! suds_t::model.loaded() )
    suds_t::model.read( param.has( "model" ) ? param.value( "model" ) : "_1" , 
			"" , "" ,
			param.has( "sig" ) && param.value( "sig" ) != "*" ? param.value( "sig") : "C4_M1"
			);
  
  suds_indiv_t self;
  self.place( edf , param , stagefile );

}

// REBASE : change epoch duration 
void proc_rebase_soap( edf_t & edf , param_t & param  )
{
  
  // uses existing EPOCH settings (which can include overlapping epochs)
  // and then predict into a new space ;  
  
  if ( ! edf.timeline.epoched() ) 
    Helper::halt( "REBASE requires that EPOCH was explicitly set beforehand" );
  
  // original epoch size (i.e. from a previous EPOCH command )
  //uint64_t e1_tp = edf.timeline.epoch_len_tp_uint64_t();
  
  // final desired epoch size (can be anything)
  double newlen = param.requires_dbl( "dur" );

  //uint64_t e2_tp = e2 * globals::tp_1sec ; 
  
  // analysis window: size must be an integer multiple of original epoch length
  // i.e. to ensure integer number of sample points (although prob not really necessary...)
  
  // double win2 = param.has( "window" ) ? param.requires_dbl( "window" ) : 0;
  // uint64_t win2_tp = win2 * globals::tp_1sec ;
  
  // double overlap2 = param.has( "overlap" ) ? param.requires_dbl( "overlap" ) : 0;
  //uint64_t overlap2_tp = overlap2 * globals::tp_1sec ;
  
  // assumption is that the set epoch length is the GCD (greatest common divisor) of the 
  // original and analysis window length
  
  // if ( win2_tp % e1_tp != 0 ) 
  //   Helper::halt( "window must be an exact multiple of current epoch length" );
  
  suds_t::set_options( param );
  
  // load model, if not already done
  if ( ! suds_t::model.loaded() )
    suds_t::model.read( param.has( "model" ) ? param.value( "model" ) : "_1" ,
                        "" , "" ,
			param.has( "sig" ) && param.value( "sig" ) != "*" ? param.value( "sig") : "C4_M1"
                        );

  suds_indiv_t self;
  self.rebase( edf , param , newlen );
  
}


// RESOAP : single observation accuracies and probabilities;
//  given a previous call to SOAP, can update observed stages
//  and model/predictions (i.e. for iterative staging)
void proc_resoap( edf_t & edf , param_t & param  )
{  
  // check that this same individual has been cached by 
  // a previous SOAP run
  if ( suds_t::cached.id != edf.id ) 
    Helper::halt( "need to SOAP w/ 'save' option before running RESOAP" );

  // need to reset only y[]
  // keep obs_stage[] and obs_stage_valid[] as is (i.e. any 'original' true staging)

  //
  // scrub all stages?
  //
  
  if ( param.has( "scrub" ) )
    {
      for (int i=0; i < suds_t::cached.y.size(); i++)
	suds_t::cached.y[i] = suds_t::str( SUDS_UNKNOWN );            
      return;
    }
  
  //
  // pick N of each epoch at random?
  //
  
  if ( param.has( "pick" ) )
    {
      int n = param.requires_int( "pick" );
      suds_t::cached.resoap_pickN( edf , n );
      suds_t::cached.resoap( edf , param.has( "verbose" ) );
      return;
    }

  //
  // else, alter a single epoch
  //
  
  // which epoch is being updated...
  int epoch = param.requires_int( "epoch" );
  // ...to which stage?
  suds_stage_t stage = suds_t::type( param.requires( "stage" ) );

  // update and refit model based on set PSC 
  suds_t::cached.resoap_alter1( edf , epoch , stage );
  suds_t::cached.resoap( edf , param.has( "verbose" ) );
  
}

// MAKE-SUDS : populate folder 'db' with trainers
void proc_make_suds( edf_t & edf , param_t & param  )
{
  // misc options
  suds_t::set_options( param );  
  
  // load model, if not already done
  if ( ! suds_t::model.loaded() ) 
    suds_t::model.read( param.requires( "model" ) );
  
  // load this individual's data, process and output text-format in 'db' folder
  suds_indiv_t trainer;
  trainer.add_trainer( edf , param );  
}


// EVAL-STAGES : given an external file (.eannot), calculate kappa and all other POPS stats
void proc_eval_stages( edf_t & edf , param_t & param )
{
#ifdef HAS_LGBM
  // one external file, vesus internal staging
  pops_indiv_t indiv( edf , param , param.requires( "file" ) );
#else
  Helper::halt( "no LGBM support compiled in" );
#endif

}

// E-CYCLE : empircal NREM cycles
void proc_ecycle(  edf_t & edf , param_t & param )
{
  // note: can relax the LGBM dependency here most likely...
#ifdef HAS_LGBM
  ecycles_t ecycles = ecycles_t( edf , param );
#else
  Helper::halt( "no LGBM support compiled in" );
#endif

}

// RUN-POPS : wrapper to run POPS
void proc_runpops( edf_t & edf , param_t & param )
{
#ifdef HAS_LGBM

  // use same syntax as lunapi 
  // sig=S1 or sig=S1,S2
  // ref=R1 or sig=R1,R2 --> do_reref
  
  // edger=Y/N  (default Y)
  // filter=Y/N (default Y)
 

  // main signal (currently, either 1 or 2) 
  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );
  if ( signals.size() == 0 ) Helper::halt( "no signals found matching " + param.value( "sig" ) );
  
  const int ns = signals.size();
  // use [x][y] format for sig specification
  
  const std::string slab = Helper::stringize( signals.signal_labels , "," );
  
  // reference channels?
  const bool do_reref = param.has( "ref" );

  // use vec instead of signal_list_t, as we may want to same ref for mult channels
  // means this will not respect aliases
  std::vector<std::string> refs;
  if ( do_reref )
    {
      refs = param.strvector( "ref" );    
      if (refs.size() != signals.size() )
	Helper::halt( "if specified, ref must match sig length" );
    }
  
  //
  // (removed for now) set LON and LOFF? (as ivars)
  //
  
  // if ( param.has( "lights-on" ) )
  //   cmd_t::ivars[ edf.id ][ "LON" ] = param.value( "lights-on" );
  
  // if ( param.has( "lights-off" ) )
  //   cmd_t::ivars[ edf.id ][ "LOFF" ] = param.value( "lights-off" );
  

  // other optional POPS args to pass in  ' args="op1=val1 op2=val2" ' 
  const std::string opt_args = param.has( "args" ) ? param.value( "args" ) : "";
    
  // ignore extant staging?
  const bool ignore_obs_staging = param.has( "ignore-obs" ) ? param.yesno( "ignore-obs" ) : false ;
  
  // required POPS library (defaults to s2)
  const std::string pops_lib = param.has( "lib" ) ? param.requires( "lib" ) : "s2";
  
  // optional path
  const std::string pops_path = param.has( "path" ) ? param.value( "path" ) : "." ; 

  // other options
  const bool do_filter = param.has( "filter" ) ? param.yesno( "filter" ) : true ;

  const bool do_edger  = param.has( "edger" ) ? param.yesno( "edger" ) : true;  


  // xsigs handles
  
  std::string allsigs = Helper::xsigs( "[" + slab + "][_F]" );
  std::string allzigs = Helper::xsigs( "[" + slab + "][_F_N]" );				      
  
  //
  // copy signal (do not alter original)
  //
  
  logger << "  ------------------------------------------------------------\n"
	 << "  making copies of original signals (appending _F)\n";

  param_t copy_param;
  copy_param.add( "sig" , slab );
  copy_param.add( "tag" , "_F" );
  proc_copy_signal( edf , copy_param );

    
  //
  // referencing?
  //

  if ( do_reref )
    {
      logger << "  ------------------------------------------------------------\n"
	     << "  re-rereferncing signals\n";

      for (int s=0; s<ns; s++)
	{
	  param_t reref_param;
	  reref_param.add( "sig" , signals.label(s) + "_F" );
	  reref_param.add( "ref" , refs[s] );
	  proc_reference( edf , reref_param );
	}
    }

  //
  // resample if needed (fixed 128 Hz)
  //

  logger << "  ------------------------------------------------------------\n"
	 << "  resampling " << allsigs << " to 128 Hz if needed\n";
  
  param_t resample_param;
  resample_param.add( "sig" , allsigs );
  resample_param.add( "sr" , "128" );
  proc_resample( edf , resample_param );
  
  //
  // filter
  //

  if ( do_filter )
    {
      logger << "  ------------------------------------------------------------\n"
	     << "  bandpass filtering signals\n";
      
      param_t filter_param;
      filter_param.add( "sig" , allsigs );
      filter_param.add( "bandpass" , "0.3,35" );
      filter_param.add( "tw" , "0.2" );
      filter_param.add( "ripple" , "0.01" );
      proc_filter( edf , filter_param );

    }

  //
  // copy signal
  //
  
  logger << "  ------------------------------------------------------------\n"
	 << "  making time-domain normalized signals\n";
  
  param_t copy2_param;
  copy2_param.add( "sig" , allsigs );
  copy2_param.add( "tag" , "_N" );
  proc_copy_signal( edf , copy2_param );
  
  //
  // normalize 
  //

  param_t norm_param;
  norm_param.add( "sig" , allzigs );
  norm_param.add( "epoch" );
  norm_param.add( "winsor" , "0.002" );
  proc_standardize( edf , norm_param );
  
  //
  // optional edger tool (on filtered signal only)
  //

  if ( do_edger )
    {
      logger << "  ------------------------------------------------------------\n"
	     << "  scanning to trim excess leading/trailing wake/artifact\n";
      
      param_t edger_param;
      edger_param.add( "sig" , allsigs );
      edger_param.add( "cache" , "ec1" );

      if ( ignore_obs_staging )
        edger_param.add( "all" );

      proc_trim( edf , edger_param );
    }
      

  //
  // run POPS
  //

  logger << "  ------------------------------------------------------------\n"
	 << "  running POPS\n";
      
  param_t pops_param;
  pops_param.add( "force-reload" );
  pops_param.add( "path" , pops_path );
  pops_param.add( "lib" , pops_lib );
  pops_param.add( "cache" , "ec1" );
  pops_param.add( "alias" , "CEN,ZEN|" + signals.label(0) + "_F," + signals.label(0) + "_F_N" );
  // equiv channels
  if ( ns > 1 )
    {
      std::string eq = "CEN,ZEN";
      
      for (int s=1; s<ns; s++)
	eq += "|" + signals.label(s) + "_F," + signals.label(s) + "_F_N";
      
      pops_param.add( "equiv" , eq );
    }

  if ( ignore_obs_staging )
    pops_param.add( "ignore-obs-staging" );

  if ( opt_args != "" )
    {
      logger << "  adding additional args to POPS: " << opt_args << "\n";
      std::vector<std::string>  tok = Helper::parse( opt_args , " \t" );
      for (int i=0; i<tok.size(); i++)
	pops_param.parse( tok[i] );
    }

  proc_pops( edf , pops_param );

  
  //
  // drop signals
  //

  logger << "  ------------------------------------------------------------\n"
	 << "  cleaning up temporary signals\n";
  
  param_t drop_param;
  drop_param.add( "sig" , allsigs + "," + allzigs );
  drop_param.add( "drop" );
  proc_drop_signals( edf , drop_param );

    
#else
  Helper::halt( "no LGBM support compiled in" );
#endif

}

// POPS : population-level staging
void proc_pops( edf_t & edf , param_t & param )
{

#ifdef HAS_LGBM

  //
  // set up parameters
  //

  pops_t pops( param );   

  //
  // force new specs? (for use w/ R/moonlight)
  //

  if ( param.has( "force-reload" ) )
    {
      pops_t::specs.init();
      pops_t::specs.init_default();     
    }

  //
  // set up features ('.' = use internal defaults)
  //
  
  std::string feature_file  = ".";
  if ( param.has( "features" ) )
    feature_file = param.value( "features" );
  else if ( pops_opt_t::pops_root != "" )
    feature_file = pops_opt_t::pops_root + ".ftr";  
  if ( feature_file != "." )
    feature_file = pops_t::update_filepath( feature_file );  
  if ( feature_file == "." ) 
    {
      Helper::halt( "POPS requires a feature file, via lib or features args" );
      return; 
    }
  
  pops_t::specs.read( feature_file );
 
  //
  // process individual (either trainer, or target)
  //
  
  pops_indiv_t indiv( edf , param );

#else
  Helper::halt( "no LGBM support compiled in" );
#endif
    
}


// SUDS : staging
void proc_suds( edf_t & edf , param_t & param )
{
  // clear cache?
  if ( param.has( "clear" ) )
    {
      suds_t::empty_banks();
      logger << "  clearing SUDS cache\n";
      return;
    }
  
  // set up global parameters (i.e. should apply to target /and/ all trainers)
  suds_t suds;
  suds_t::set_options( param );

  // load model, if not already done                                                                                                     
  if ( ! suds_t::model.loaded() )
    suds_t::model.read( param.requires( "model" ) );

  // load trainers, if not already done
  //
  // bank() and wbank() can share the same individuals 
  // (but nore, they are loaded twice, currently... need to fix)
  //
  // The feature matrix X is only loaded for wdb members
  
  if ( param.has( "wdb" ) )
    {
      // load as separate files (i.e. duplicate) 
      suds.attach_db( param.requires( "db" ) , true , false );
      suds.attach_db( param.value( "wdb" ) , false , true );
    }
  else if ( param.has( "db" ) )
    {
      suds.attach_db( param.value( "db" ) , true , false );
    }
  else if ( param.has( "lib" ) ) 
    {
      suds.attach_lib( param.value( "lib" ) );
    }
  else
    Helper::halt( "no library attached" );
	   
  // do actual scoring  
  suds.score( edf , param );
  
}




// ZR : Z-ratio 
void proc_zratio( edf_t & edf , param_t & param )
{    
  std::string signal = param.requires( "sig" );

  staging_t staging;
  staging.zratio.calc( edf , signal );
}


// CORRECT : use regression or PCA to correct artifacts 

void proc_correct( edf_t & edf , param_t & param )
{
  dsptools::artifact_correction( edf , param );
}

// COMBINE : make a new signal
void proc_combine( edf_t & edf , param_t & param )
{
  edf.combine( param );
}

// TRIM : empircal heuristic to set lights off/on times based on signals
void proc_trim( edf_t & edf , param_t & param )
{
  dsptools::trim_lights( edf, param );
}

// ARTIFACTS : artifact rejection using Buckelmueller et al. 
void proc_artifacts( edf_t & edf , param_t & param )	  
{
  std::string signal = param.requires( "sig" );
  annot_t * a = buckelmuller_artifact_detection( edf , param , signal );  
}

// MOVING-AVERAGE
void proc_moving_average( edf_t & edf , param_t & param )
{
  dsptools::movavg( edf , param );  
}

// FILTER : general FIR
void proc_filter( edf_t & edf , param_t & param )	  
{

  if ( param.has( "butterworth" ) || param.has( "chebyshev" ) )
    dsptools::apply_iir( edf , param );
  else
    dsptools::apply_fir( edf , param );
}


// -fir  from the command line
void proc_filter_design_cmdline()
{
  
  // expect parameters on STDIN
  
  param_t param;
  while ( ! std::cin.eof() )
    {
      std::string x;
      std::cin >> x;      
      if ( std::cin.eof() ) break;
      if ( x == "" ) continue;
      param.parse( x ); 
    }

  dsptools::design_fir( param );
  
}

// TV   total variation 1D denoising
void proc_tv_denoise( edf_t & edf , param_t & param )
{
  dsptools::tv( edf , param );
}

// OTSU   automatic binary thresholding
void proc_otsu( edf_t & edf , param_t & param )
{
  dsptools::otsu( edf , param );
}

// CWT 
void proc_cwt( edf_t & edf , param_t & param )
{
  dsptools::cwt( edf , param );
}

// HILBERT 
void proc_hilbert( edf_t & edf , param_t & param )
{
  dsptools::hilbert( edf , param );
}

// SYNC
void proc_sync( edf_t & edf , param_t & param )
{
  dsptools::sync( edf , param );
}

// TSYNC
void proc_tsync( edf_t & edf , param_t & param )
{
  dsptools::tsync( edf , param );
}

// XCORR - better TSYNC for basic XCORR
void proc_xcorr( edf_t & edf , param_t & param )
{
  dsptools::xcorr( edf , param );
}

// -cwt  from the command line
void proc_cwt_design_cmdline()
{
  
  // expect parameters on STDIN
  
  param_t param;
  while ( ! std::cin.eof() )
    {
      std::string x;
      std::cin >> x;      
      if ( std::cin.eof() ) break;
      if ( x == "" ) continue;
      param.parse( x ); 
    }
  
  dsptools::design_cwt( param );
  
}


// -copy-suds-db  from the command line
void proc_copy_suds_cmdline()
{

  // this takes a SINGLE text-format file,
  // (which may contain multiple individuals)
  // and writes out a SINGLE binary file

  // we now have the following functions only
  //  MAKE-SUDS   : write a library ( text format, one file per trainer only, written to a folder)
  //  cat         : merge library (multiple individuals to single library)
  //  --copy-suds : reformat text->binary single 
  //  SUDS        : read a single binary format file
  
  // expect parameters on STDIN
  param_t param;
  while ( ! std::cin.eof() )
    {
      std::string x;
      std::cin >> x;      
      if ( std::cin.eof() ) break;
      if ( x == "" ) continue;
      param.parse( x ); 
    }  

  std::string f1 = param.requires( "from" );
  std::string f2 = param.requires( "to" );
  suds_t::text2binary( f1 , f2 , param.has( "with-features" ) ) ;
  
}

// --combine-suds  from the command line
void proc_combine_suds_cmdline()
{

  // this takes a SINGLE text-format file,
  // (which may contain multiple individuals)
  //  1) extracts the features and stages 
  //  2) creates a single 'collated' individual (based on the existing features)
  //  3) runs SVD as usual
  //  4) writes the output as a single text-format file
  //   (i.e.  which can then be transformed w/ --copy-suds)

  // we now have the following functions only
  //  MAKE-SUDS   : write a library ( text format, one file per trainer only, written to a folder)
  //  cat         : merge library (multiple individuals to single library)
  //  --combine-suds : merge feature matrics across individuals and recompute SVD
  //  --copy-suds : reformat text->binary single  
  //  SUDS        : read a single binary format file
  
  // expect parameters on STDIN
  param_t param;
  while ( ! std::cin.eof() )
    {
      std::string x;
      std::cin >> x;      
      if ( std::cin.eof() ) break;
      if ( x == "" ) continue;
      param.parse( x ); 
    }  

  suds_t::combine_trainers( param );
  
}

// FILTER-DESIGN : general FIR design
void proc_filter_design( edf_t & edf , param_t & param )	  
{
  dsptools::design_fir( param );
}

// CWT-DESIGN : CWT design
void proc_cwt_design( edf_t & edf , param_t & param )	  
{
  dsptools::design_cwt( param );
}

// ZOH : special case of upsampling
void proc_zoh( edf_t & edf , param_t & param )
{
  dsptools::resample_channel_zoh( edf, param );
}

// RESAMPLE : generic sample-rate conversion 
void proc_resample( edf_t & edf , param_t & param ) 
{
  dsptools::resample_channel( edf, param );
}

// FIX-SAMPLING : cublic spline to fix bad EDF record/SR combos
void proc_fix_sampling( edf_t & edf , param_t & param )
{
  dsptools::fix_sampling( edf , param );
}

// MS: microstate analysis
void proc_microstates( edf_t & edf , param_t & param )
{
  dsptools::microstates( edf , param );
}

// ASYMM
void proc_asymm( edf_t & edf  , param_t & param )
{
  lat_t lat( edf , param );
}

// PERI
void proc_peri( edf_t & edf  , param_t & param )
{
  // generic peri-event epoch-based analyses
  dsptools::peri( edf , param );
}

// TLOCK
void proc_tlock( edf_t & edf  , param_t & param )
{
  // get mean time-locked value of one signal against a set of annotations (time-points)
  dsptools::tlock( edf , param );
}

// TCLST
void proc_tclst( edf_t & edf  , param_t & param )
{
  // designed for SW analysis, given a set of points, cluster the channels/signals around them  
  dsptools::tclst( edf , param );
}

// KMEANS - apply K-means to time-series
//  - assumes will be applied in low sample count context, e.g.
//    on 1 Hz signals, etc 
void proc_kmeans( edf_t & edf , param_t & param )
{
  //  dsptools::kmeans( edf , param );
}

// PEAKS
void proc_peaks( edf_t & edf , param_t & param )
{
  dsptools::peaks( edf , param );
}

// Z-PEAKS
void proc_zpeaks( edf_t & edf , param_t & param )
{
  dsptools::zpeaks( edf , param );
}

// SEDF : make a summarize EDF 
void proc_sedf( edf_t & edf , param_t & param )
{
  sedf_t sedf( edf , param );  
}


// PREDICT : given cached values, make a model-based prediction
void proc_predict( edf_t & edf , param_t & param )
{
  prediction_t m( edf , param );
}

// PSC : either build PSC (from multiple results) or fit to an EDF
void proc_psc( edf_t & edf , param_t & param )
{

  // optionally, we can clear the attached W/V store, e.g. if 
  // applying the same cached metrics to multiple projections: 
  
  // PSI sig=... cache-metrics=c1
  // TAG P/1
  // PSC cache=p1 proj=proj1.txt 
  // PSC clear
  // TAG P/2
  // PSC cache=p1 proj=proj2.txt 

  if ( param.has( "clear" ) )
    psc_t::clear_proj();
  
  // note:  PSC contruct() is called w/out an EDF from the command line
  //  luna --psc <args>

  psc_t psc;
  
  // if already populated, this returns so we can call multiple times
  // expects 'proj=' variable to point to output of 'proj=' from PSC 
  psc.attach( param );
  
  // project attached solution
  psc.project( edf , param );

}

// DYNAM : take arbitrary inputs (files/ signals)
//         and run qdynam_t
void proc_qdynam( edf_t & edf , param_t & param )
{
  // assumes a) has HYPNO done
  //         b) inputs will match current epoching
  
  dsptools::qdynam( edf , param );
  
}


// PSD : calculate PSD via Welch
void proc_psd( edf_t & edf , param_t & param )	  
{  
  std::string signal = param.requires( "sig" );
  annot_t * power = spectral_power( edf , signal , param );  
}

// FFT : caclulate basic FFT
void proc_fft( edf_t & edf , param_t & param )
{
  dsptools::fft( edf , param );
}

// MTM : calculate MTM 
void proc_mtm( edf_t & edf , param_t & param )	  
{  
  mtm::wrapper( edf , param );
}

// 1FNORM : normalization of signals for the 1/f trend
void proc_1overf_norm( edf_t & edf , param_t & param )	  
{  
  dsptools::norm_1overf( edf,  param ) ;
}

// IRASA : calculate IRASA
void proc_irasa( edf_t & edf , param_t & param )
{  
  irasa_wrapper( edf , param );
}

// FIP: frequency/interval plot
void proc_fiplot( edf_t & edf , param_t & param )	  
{  
  fiplot_wrapper( edf , param );
}

// TAG : analysis tag
void proc_tag( param_t & param )
{
  // either TAG tag=lvl/fac
  // or just TAG lvl/fac 

  if ( ! param.single() ) Helper::halt( "TAG requires a single argument" );

  if ( param.has( "tag" ) )
    set_tag( param.value( "tag" ) );
  else 
    set_tag( param.single_value() );
  
}

void set_tag( const std::string & t ) 
{
  globals::current_tag = t ; 

  if ( t != "." ) 
    logger << "  setting analysis tag to [" << globals::current_tag << "]\n";

  if ( t == "." ) writer.tag( "." , "." );
  else
    {
      std::vector<std::string> tok = Helper::parse( globals::current_tag , "/" );
      if ( tok.size() != 2 ) Helper::halt( "TAG format should be factor/level" );

      // check that it is not a known tage
      std::string fac = Helper::toupper( tok[0] );

      if ( fac == globals::freq_strat ||
	   fac == globals::signal_strat ||
	   fac == globals::stage_strat  ||
	   fac == globals::cycle_strat  ||
	   fac == globals::band_strat   ||
	   fac == globals::annot_strat ||
	   fac == globals::annot_instance_strat ||
	   fac == globals::annot_meta_strat ||
	   fac == globals::count_strat  ||
	   fac == globals::epoch_strat  ||
	   fac == globals::time_strat  ||
	   fac == globals::sample_strat ||
	   fac == globals::cluster_strat ||
	   fac == "TH" || fac == "MSEC" || fac == "SP" )
	Helper::halt( "cannot use " + tok[0] + " as a TAG factor, matches an internal label" );
            
      // set
      writer.tag( tok[1] , tok[0] );
    }
}

// ANON : anonymize EDF 
void proc_anon( edf_t & edf , param_t & param )
{

  // either
  //  1.  Set EDF ID/date to NULL (respecting silly EDF+ conventions)
  //  2.  Set EDF ID to (sample list) ID  ('insert-id')
  //  3.  Change EDF /and/ (sample-list) ID to root_000000N  ('root')

  const std::string anon_id = edf.header.edfplus ? "X X X X" : ".";
  const std::string rec_info = edf.header.edfplus ? "Startdate X X X X" : ".";

  if ( param.has ( "insert-id" ) ) 
    {
      logger << " setting ID to " << edf.id << " and start date to '01.01.85' for " 
	     << edf.filename << "\n";
      
      edf.header.patient_id = edf.header.edfplus ? edf.id + " X X X" : edf.id;
      
    }
  else if ( param.has( "root" ) )
    {

      ++globals::anon_idroot_cnt;

      // root_N
      std::string newid = param.value( "root" ) + "_" + Helper::int2str( globals::anon_idroot_cnt );
      
      edf.header.patient_id = edf.header.edfplus ? newid + " X X X" : newid;
      
      edf.id = newid;

      logger << " setting ID and EDF ID to " << newid << "\n";

    }
  else
    {
      logger << " setting ID and start date to null ('" << anon_id << "' and '01.01.85') for " 
	     << edf.filename << "\n";

      edf.header.patient_id = anon_id;

    }

  // recording info field
  edf.header.recording_info = rec_info;

  // 'clipping' date for EDF  
  edf.header.startdate = "01.01.85";
}


// DUMP : dump all data
void proc_dump( edf_t & edf , param_t & param )	  
{
  std::string signal = param.requires( "sig" );  
  edf.data_dumper( signal , param );	  
}


// INSERT
void proc_insert( edf_t & edf , param_t & param )
{
  edf_inserter_t inserter( edf , param );
}

// ALIGN-EPOCHS
void proc_align_epochs( edf_t & edf , param_t & param )
{
  align_epochs_t align( edf , param );  
}

// ALIGN-ANNOTS: given an ALIGN-EPOCHS solution, change annotation timings
void proc_align_annots( edf_t & edf , param_t & param )
{
  align_annots_t align( edf , param );  
}

// ALIGN
void proc_align( edf_t & edf , param_t & param )
{

  // requires one or more annotation codes;
  // find only annotations that completely span the signal (i.e. taking discontinuous EDFs into account)
  // remove samples that aren't spanned
  // align remaining samples to be in whole EDF records, 
  // and rewrite the EDF, i.e. like 

  // only pull out samples that overlap with these annotations, where a 

  if ( ! param.has( "align" ) ) Helper::halt( "no 'align' annotations specified" );
  
  std::vector<std::string> a = param.strvector( "align" );

  logger << "  realigning EDF based on annotation list: " << param.value( "align" ) << "\n";

  bool okay = edf.align( a );
  
  if ( ! okay ) 
    {
      logger << "  problem in creating the aligned EDF, bailing...\n"
	     << "  (check there are 1+ valid channels)\n";      
      return;
    }
  
  logger << "  now WRITE'ing realigned EDF (and annotations if 'annot-out' set) to disk\n"
	 << "  note:  this will will set the 'problem' flag to skip to next EDF\n";
  
  proc_write( edf , param );
  
  if ( param.has( "annot-out" ) )
    edf.annotations->write( param.requires( "annot-out" ) , param , edf );
  
  // force this to be the last command executed for this EDF, i.e. no it is saved
  globals::problem = true;

}


// EPOCH DUMP 
void proc_epoch_dump( edf_t & edf , param_t & param )
{
  // REDUNDANT ; command not documented
  std::set<std::string> * annots = NULL;
  if ( param.has( "annot" ) )
    {
      annots = new std::set<std::string>;
      *annots = param.strset( "annot" ); // parse comma-delim list
    }

  edf.data_epoch_dumper( param , annots );
}


// MATRIX 
void proc_epoch_matrix( edf_t & edf , param_t & param )
{  
  edf.epoch_matrix_dumper( param );
}


// HEAD
void proc_head_matrix( edf_t & edf , param_t & param )
{
  edf.head_matrix_dumper( param );
}


// INTERVALS : raw signal data from an interval list
void proc_intervals( param_t & param , const std::string & data )	  
{  

  std::string ints = param.requires( "intervals" );
  // e.g.: INTERVAL edfs=../scratch/edf.list 
  // intervals from STDIN
  dump_intervals( ints , data );
}


// PSI : phase slope index
void proc_psi( edf_t & edf , param_t & param )
{
  dsptools::psi_wrapper( edf , param );
}

// COVAR : covariance between two signals (not implemented)
void proc_covar( edf_t & edf , param_t & param )
{  
  std::string signals1 = param.requires( "sig1" );
  std::string signals2 = param.requires( "sig2" );
  edf.covar(signals1,signals2);
}

// AROUSALS : detect micro-arousals during sleep 
void proc_arousals( edf_t & edf , param_t & param )
{
  arousals_t arousals( edf , param );  
}

// SPINDLES : spindle detection using CWT or bandpass/RMS
void proc_spindles( edf_t & edf , param_t & param )
{	

  // default to wavelet analysis
  std::string method = param.has( "method" ) ? param.value( "method" ) : "wavelet" ; 
  
  annot_t * a = NULL;
  
  if ( method == "wavelet" ) a = spindle_wavelet( edf , param );
  else Helper::halt( "SPINDLE method not recognized; should be 'bandpass' or 'wavelet'" );

}

// COUPL : spindle/SO couplig
void proc_coupling( edf_t & edf , param_t & param )
{
  // requires cached SPINDLES and SO results
  spindle_so_coupling( edf , param );
}

// PCOUPL : generic phase coupling
void proc_generic_coupling( edf_t & edf , param_t & param )
{
  // requires cached SPINDLES and SO results
  dsptools::phase_coupling( edf , param );
}

// RIPPLES : ripple detection
void proc_ripples( edf_t & edf , param_t & param )
{
  dsptools::ripple_wrapper( edf , param );
} 

// POL : polarity check for EEG N2/N3 
void proc_polarity( edf_t & edf , param_t & param )
{	
  dsptools::polarity( edf , param );
}

// QC : generic PSG QC scan
void proc_qc( edf_t & edf , param_t & param )
{
  dsptools::qc_t qc( edf , param );
}

// REMS : detect REMS via simple heuristic
void proc_rems( edf_t & edf , param_t & param )
{
  dsptools::rems( edf , param );
}

// SO : detect slow waves, do time-locked FFT on rest of signal
void proc_slowwaves( edf_t & edf , param_t & param )
{	
  slow_waves_t sw( edf , param );  
}


// EDF-MINUS : convert from EDF+D to EDF
//   adding padding (zeros for annots)
//   ajusting annotations
//   and adding in a new "gap" annot
void proc_edf_minus( edf_t & edf , param_t & param )
{
  edf.edf_minus( param );
}

// SET-TIMESTAMPS
void proc_set_timestamps( edf_t & edf , param_t & param )
{
  edf.set_timestamps( param );
}

// EDF : convert from EDF+D to EDF or EDF+C
//               or EDF+C to EDF
void proc_force_edf( edf_t & edf , param_t & param )
{

  Helper::halt( "DROP-TIMESTAMPS is not available" );
  
  if ( ! edf.header.edfplus ) 
    {
      logger << "  already a standard EDF, nothing to do\n";
      return;
    }

  if ( edf.header.continuous )
    {
      logger << "  converting from EDF+C to standard EDF\n";
      edf.set_edf();      
      edf.reset_start_time();      
      edf.timeline.set_epoch( globals::default_epoch_len , globals::default_epoch_len );
      edf.timeline.re_init_timeline();
      edf.restructure( true );
      return;
    }
  
  // if here, it is nominally EDF+D
  // working point here.

  if ( ! edf.is_actually_discontinuous() )
    {
      logger << "  converting from EDF+D that is actually continuous, to standard EDF\n";
      edf.set_edf();
      edf.reset_start_time();      
      edf.timeline.re_init_timeline();
      edf.timeline.set_epoch( globals::default_epoch_len , globals::default_epoch_len );
      return;
    }
  
  //
  // actual EDF+D -> EDF conversion (drop time-line)
  //

  logger << "  forcing EDF+D to standard EDF: will lose discontinuity/time information\n";  
  edf.set_edf();
  
  // set start time to NULL
  logger << "  setting EDF starttime to null (00.00.00)\n";
  edf.header.starttime = "00.00.00";
  
  edf.timeline.set_epoch( globals::default_epoch_len , globals::default_epoch_len );
  edf.timeline.re_init_timeline();
  edf.restructure( true );
  return;
  
}


// WRITE : write a new EDF or EDFZ to disk
// optionally, writes annotation files,
// (although user responsible for not altering time structure of EDF
// i.e. we do not change the time encoding of the annotation file, which
// are always anchored to the original)
void proc_write( edf_t & edf , param_t & param )
{
  
  // write a .edfz and .edfz.idx
  const bool edfz = param.yesno( "edfz" );
  
  // add 'tag' to new EDF
  std::string filename = edf.filename;

  if ( Helper::file_extension( filename, "edf" ) || 
       Helper::file_extension( filename, "EDF" ) ) 
    filename = filename.substr(0 , filename.size() - 4 );
  
  if ( Helper::file_extension( filename, "edfz" ) || 
       Helper::file_extension( filename, "EDFZ" ) ) 
    filename = filename.substr(0 , filename.size() - 5 );

  if ( Helper::file_extension( filename, "edf.gz" ) || 
       Helper::file_extension( filename, "EDF.GZ" ) ) 
    filename = filename.substr(0 , filename.size() - 7 );

  // make edf-tag optional
  if ( param.has( "edf" ) )
    filename = param.requires( "edf" ) + ".edf" ;
  else if ( param.has( "edf-tag" ) ) 
    filename += "-" + param.requires( "edf-tag" ) + ".edf";
  else
    {
      // but here, must specify a folder explicitly in that case...
      if ( ! param.has( "edf-dir" ) )
	Helper::halt( "if not adding edf-tag, must explicitly specify edf-dir" );
      
      filename += ".edf";
    }

  // set .edf.gz as main EDFZ file extension
  if ( edfz ) filename += ".gz";
  
  //
  // optionally, allow directory change
  //

  if ( param.has( "edf-dir" ) )
    {
      
      std::string outdir = param.value("edf-dir");
      
      if ( outdir[ outdir.size() - 1 ] != globals::folder_delimiter ) 
	outdir += globals::folder_delimiter;

      //Helper::halt("edf-dir value must end in '" + std::string(1,globals::folder_delimiter) + "' to specify a folder" );

      int p=filename.size()-1;
      int v = 0;
      for (int j=p;j>=0;j--)
	{
	  if ( filename[j] == globals::folder_delimiter ) { v=j+1; break; }
	}
      filename = outdir + filename.substr( v );            
      
      // create folder if it does not exist 
      // -p is (usually?) not needed for Windows

      std::string syscmd = globals::mkdir_command + " " + param.value( "edf-dir" );

      int retval = system( syscmd.c_str() );
      
    }

  //
  // Sample list
  //

  if ( param.has("sample-list") )
    {	  
      std::string file = param.value("sample-list");
      
      // also append annotation files
      bool append_annots = param.has( "with-annots" );

      // open/append
      logger << "  appending " << filename
	     << " to sample-list " << file
	     << ( append_annots ? " (with annotations)" : " (dropping any annotations)" ) << "\n";
      
      std::ofstream FL( file.c_str() , std::ios_base::app );
      FL << edf.id << "\t"
	 << filename;
      if ( append_annots )
	for (int i=0;i<edf.annot_files.size();i++) FL << "\t" << edf.annot_files[i];
      FL << "\n";
      
      FL.close();
    }

  //
  // prep EDF for writing then write to disk
  //

  // arbitrary, but need epochs if not set it seems
  if ( !edf.timeline.epoched() ) 
    edf.timeline.set_epoch( 30 , 30 );

  // if a mask has been set, this will restructure the mask
  edf.restructure(); 
  

  //
  // Force as EDF (i.e. even if restructured), and set starttime = 0;
  //

  int write_as_edf = param.has( "force-edf" ) ? 1 : 0 ;

  if ( param.has( "null-starttime" ) )
    {
      if ( ! write_as_edf )
	Helper::halt( "null-starttime option can only be specified with force-edf" );
      write_as_edf = 2;
    }  
  
  
  //
  // Do not write 'quasi-discontinuous' (i.e. single-segment EDF+D) as EDF+D
  // Rather, write as EDF+C
  //

  bool always_EDFD = param.has( "EDF+D" );

  //
  // Specify order and number of channels in the new EDF?
  //
  
  std::vector<int> channels;
  
  const bool set_chorder = param.has( "channels" );
  
  if ( set_chorder )
    {
      std::vector<std::string> str = param.strvector( "channels" );
      std::set<int> cuniq;
      // check each exists
      for (int s=0; s<str.size(); s++)
	{
	  if ( ! edf.header.has_signal( str[s] ) ) 
	    Helper::halt( "could not find requested channel " + str[s] );
	  const int slot = edf.header.signal( str[s] );
	  channels.push_back( slot );
	  cuniq.insert( slot );
	}

      if ( cuniq.size() < channels.size() )
	logger << "  exporting " << cuniq.size()
	       << " unique signals ("<< channels.size()
	       << " total) from " << edf.header.ns << " originals\n";
      else
	logger << "  exporting " << channels.size() << " signals from "	       
	       << edf.header.ns << " originals\n";
    }
  
  // Save data (write_as_edf flag forces starttime to 00.00.00 if
  // really EDF+D)
  
  bool saved = edf.write( filename , edfz , write_as_edf , always_EDFD , set_chorder ? &channels : NULL );
  
  if ( ! saved ) 
    Helper::halt( "problem trying to save " + filename );

}


// EPOCH : set epochs 
void proc_epoch( edf_t & edf , param_t & param )
{

  // just dump, do not alter?
  
  if ( param.has( "dump" ) || param.has( "table" ) )
    {
      
      const bool show_masked = param.has( "masked" ) ;
      
      if ( edf.timeline.epoched() )
	{
	  if ( show_masked ) 
	    logger << "  outputting epoch table for "
		   << edf.timeline.num_total_epochs()
		   << " masked & unmasked epochs\n";
	  else
	    logger << "  outputting epoch table for "
		   << edf.timeline.num_epochs()
		   << " unmasked epochs\n";
	  
	  edf.timeline.output_epoch_info( true , show_masked );
	}
      else
	logger << "  no epochs set, not dumping any information\n";
      return;
    }

  //
  // annotation-based (generic) epochs?
  //
  
  if ( param.has( "annot" ) )
    {
      int ne = edf.timeline.calc_epochs_generic_from_annots( param );
      
      logger << "  set " << ne
	     << " generic epochs, based on annotations [ " << param.value( "annot" ) << " ]";
      if ( param.has( "else" ) )
	logger << " and [ " << param.value( "else" ) << " ]";
      logger << "\n";

      edf.timeline.output_epoch_info( param.has( "verbose" ) );
      
      const bool opt_req   = param.has( "require" );
      if ( opt_req ) 
	{
	  int r = param.requires_int( "require" );
	  if ( edf.timeline.num_epochs() < r )
	    {
	      logger << " ** warning for "
		     << edf.id << " when setting EPOCH: "
		     << "required=" << r << "\t"
		     << "but observed=" << edf.timeline.num_epochs()  << "\n";
	      globals::problem = true;
	    }	  
	}
      // all done for generic epoch detection
      return;
    }

  //
  // otherwise, define standard epochs 
  //  

  const bool opt_clear = param.has( "clear" );
  const bool opt_req   = param.has( "require" );
  const bool opt_len   = param.has( "len" ) || param.has( "dur" ) || param.has( "epoch" ) || param.has( "inc" );
  const bool opt_align = param.has( "offset" ) || param.has( "align" );
  const bool span_gaps = param.has( "splice" ) || param.has( "splice-gaps" );
    
  //
  // just check we have the required number of epochs, but otherwise do not touch/re-epoch?
  // i.e. this keeps the epoch annotations intact
  //
  
  if ( opt_req && ! ( opt_clear || opt_len || opt_align ) )
    {
      int r = param.requires_int( "require" );
      if ( edf.timeline.num_epochs() < r ) 
	{	  
	  logger << " ** warning for "  
		 << edf.id << " when setting EPOCH: "
		 << "required=" << r << "\t"
		 << "but observed=" << edf.timeline.num_epochs()  << "\n";
	  globals::problem = true;
	}
      return;
    }
    
  // unepoch?
  if ( opt_clear )
    {
      logger << "  clearing all epochs: signals are now unepoched\n";
      edf.timeline.unepoch();
      return;
    }
  
  double dur = 0 , inc = 0;
  
  // default = 30 seconds, non-overlapping
  if ( ! opt_len )
    {
      dur = 30; inc = 30;
    }
  else
    {
      
      if ( param.has( "epoch" ) )
	{
	  std:: string p = param.requires( "epoch" );
	  std::vector<std::string> tok = Helper::parse( p , "," );
	  if ( tok.size() > 2 || tok.size() < 1 ) Helper::halt( "expcting epoch=length{,increment}" );
	  
	  if ( ! Helper::str2dbl( tok[0] , &dur ) ) Helper::halt( "invalid epoch length" );
	  if ( tok.size() == 2 ) 
	    {
	      if ( ! Helper::str2dbl( tok[1] , &inc ) ) 
		Helper::halt( "invalid epoch increment" );
	    }
	  else inc = dur;
	}
      
      else if ( param.has( "len" ) ) 
	{
	  dur = param.requires_dbl( "len" );	  
	  if ( param.has( "inc" ) ) 
	    inc = param.requires_dbl( "inc" );
	  else 
	    inc = dur;
	}

      else if ( param.has( "dur" ) ) 
	{
	  dur = param.requires_dbl( "dur" );	  
	  if ( param.has( "inc" ) ) 
	    inc = param.requires_dbl( "inc" );
	  else 
	    inc = dur;
	}

    }

  if ( param.has( "inc" ) )
    {
      inc = param.requires_dbl( "inc" );
    }


  //
  // Epochs start at 0, or something else?
  //
  
  uint64_t offset = 0LLU;
  
  if ( param.has( "offset" ) )
    {
      std::string ostr = param.value( "offset" );
      std::vector<std::string> tok = Helper::parse( ostr  , ":" );
      
      // hh:mm, hh:mm:ss or dd:hh:mm:ss
      // (can be hh:mm:ss.ssss)                                         
      bool is_hms = tok.size() == 2 || tok.size() == 3 || tok.size() == 4;

      if ( is_hms )
	{

	  clocktime_t starttime( edf.header.starttime );
	  
	  if ( ! starttime.valid )
	    Helper::halt( "specifying offset=hh:mm:ss clocktime start, but no valid EDF header starttime" );
	  
	  clocktime_t otime( ostr );
	  
	  // 1: EDF start comes before OFFSET start (required)
	  // 2: OFFSET comes before EDF start --> flag error	  
	  
	  int earlier = clocktime_t::earlier( starttime , otime );
	  
	  if ( earlier == 2 )
	    Helper::halt( "cannot specify an EPOCH offset earlier than EDF start" );
	  else
	    offset = globals::tp_1sec * clocktime_t::difference_seconds( starttime , otime ) ;
	  
	}
      else
	{
	  // arg value is in seconds
	  double o_sec = param.requires_dbl( "offset" ) ;
	  if ( o_sec < 0 ) Helper::halt( "offset must be non-negative" );
	  offset = globals::tp_1sec * o_sec;
	}
      
    }


  //
  // Align w/ first instance of some annotation?
  //
  
  std::vector<std::string> align_annots;
  std::string align_str = "";

  if ( param.has( "align" ) )
    {
      if ( param.has( "offset" ) ) 
	Helper::halt( "cannot specify both offset and align" );
      
      // swap in a default?
      if ( param.empty( "align" ) ) // i.e. no arg to align
	{
	  align_str = "N1,N2,N3,R,W,?,L,U,M";
	  align_annots = Helper::parse( align_str , "," );
	}
      else
	{
	  align_str = param.value( "align" );
	  align_annots = param.strvector( "align" );
	}

      // for EDF+D, this vector is passed to timeline
      
      // find the first of these annotations
      offset = edf.annotations->first( align_annots );
    }
  

  // if already EPOCH'ed for a different record size, or increment,
  // then we should first remove all epochs; this will remove the
  // EPOCH-to-record mapping too
  
  if ( edf.timeline.epoched() 
       && ( ( ! Helper::similar( edf.timeline.epoch_length() , dur ) )
	    || ( ! Helper::similar( edf.timeline.epoch_inc() , inc ) )
	    || ( ! Helper::similar( edf.timeline.epoch_offset() , globals::tp_duration * offset ) )
	    || ( edf.timeline.align_string() != align_str ) ) )
    {
      logger << " epoch definitions have changed: original epoch mappings will be lost\n";
      edf.timeline.unepoch();
    }

  //
  // gap spanning?
  //

  edf.timeline.set_epochs_to_span_gaps ( span_gaps );
  
  //
  // basic log info
  //
  
  
  int ne = edf.timeline.set_epoch( dur , inc , offset , 
				   align_str , align_annots.size() == 0 ? NULL : &align_annots );  
				   
  

  //
  // minimal output to stdout
  //

  if ( param.has( "min" ) )
    {
      std::cout << ne << "\n";
      return;
    }
  
  logger << "  set epochs, length " << dur 
	 << " (step " << inc 
	 << ", offset " << offset * globals::tp_duration  
	 <<  "), " << ne << " epochs\n";
  if ( span_gaps ) logger << "  epochs defined to splice out any gaps\n";
  
  //
  // write more verbose information to db
  //
  
  edf.timeline.output_epoch_info( param.has( "verbose" ) );

 
  //
  // any constraints on the min. number of epochs required? 
  //

  if ( param.has( "require" ) )
    {
      int r = param.requires_int( "require" );
      if ( ne < r ) 
	{	  
	  logger << " ** warning for "  
		 << edf.id << " when setting EPOCH: "
		 << "required=" << r << "\t"
		 << "but observed=" << ne << "\n";
	  globals::empty = true;
	}
    }
}


// FILE-MASK : apply a mask from a file
void proc_file_mask( edf_t & edf , param_t & param )
{ 
  std::string f = "";
  bool exclude = true;

  if      ( param.has("include") ) { f = param.requires("include"); exclude = false; }
  else if ( param.has("exclude") ) f = param.requires("exclude"); 
  else Helper::halt( "need either include or exclude for MASK-FILE");
  
  if ( param.has( "intervals" ) )
    edf.timeline.load_interval_list_mask( f , exclude );
  else
    edf.timeline.load_mask( f , exclude );
}


// EPOCH-MASK  : based on epoch-annotations, apply mask
void proc_epoch_mask( edf_t & edf , param_t & param )
{
  Helper::halt( "EPOCH-MASK command is redundant" );
}


// FREEZE : make a copy of the current (internal) EDF 
void proc_freeze( edf_t & edf , param_t & param )
{
  
  if ( ! param.single() )
    Helper::halt( "FREEZE requires a single argument" );

  // take either FREEZE tag=name or  FREEZE name
  const std::string freeze_name = param.has( "tag" ) ? param.value( "tag" ) : param.single_value() ;

  if ( freeze_name == "remove" ) Helper::halt( "cannot use 'remove' as a freeze name" );  
  
  freezer.freeze( freeze_name , edf );

}

// CLEAN : clean out entire freezer
void proc_clean_freezer( edf_t & edf , param_t & param )
{
  logger << "  cleaning the freezer...\n";
  freezer.clean( &edf );
}

// THAW : bring back and replace current EDF
void proc_thaw( edf_t & edf , param_t & param )
{
  
  const bool preserve_cache = param.has( "preserve-cache" ) ? param.yesno( "preserve-cache" ) : false ;
  
  const bool remove = param.has( "remove" ) ? param.yesno( "remove" ) : false ;
  
  const bool strictmode = param.has( "strict" ) ? param.yesno( "strict" ) : false ;
  
  bool retval;
  
  if ( remove )
    retval = freezer.thaw( param.requires( "tag" ) , &edf , remove , preserve_cache );  
  else
    {
      // can allow single arg context here (if not also using 'remove' or 'replace-cache' )
      const std::string freeze_name = param.has( "tag" ) ? param.value( "tag" ) : param.single_value() ;
      retval = freezer.thaw( freeze_name , &edf , false , preserve_cache );
    }
  
  if ( strictmode && ! retval ) 
    Helper::halt( "could not thaw requsted data freeze; under strict-mode, halting ");

}

// EPOCH-ANNOT : directly apply epoch-level annotations from the command line
// with recodes
void proc_file_annot( edf_t & edf , param_t & param )
{ 
  
  std::string f = param.requires( "file" );

  std::vector<std::string> a;
  
  std::map<std::string,std::string> recodes;
  if ( param.has( "recode" ) )
    {
      std::vector<std::string> tok = Helper::quoted_parse( param.value( "recode" ) , "," );

      for (int i=0;i<tok.size();i++)
	{
	  std::vector<std::string> tok2 = Helper::quoted_parse( tok[i] , "=" );
	  if ( tok2.size() == 2 )
	    {
	      logger << "  remapping from " << tok2[0] << " to " << tok2[1] << "\n";
	      recodes[ Helper::unquote( tok2[0] ) ] = Helper::unquote( tok2[1] );
	    }
	  else
	    Helper::halt( "bad format for " + tok[i] );
	}
    }
  
  if ( ! Helper::fileExists( f ) ) Helper::halt( "could not find " + f );
  
  std::set<std::string> amap;

  std::ifstream IN1( f.c_str() , std::ios::in );
  while ( ! IN1.eof() )
    {
      std::string x;
      Helper::safe_getline( IN1 , x );
      if ( IN1.eof() ) break;
      if ( x == "" ) continue;
      if ( recodes.find(x) != recodes.end() ) 
	{	  
	  x = recodes[x];
	}
      a.push_back( x );      
      amap.insert( x );
    }
  IN1.close();
  
  logger << " mapping " << amap.size() << " distinct epoch-annotations (" << a.size() << " in total) from " << f << "\n";


  //
  // Check correct number of epochs
  //
  
  if ( a.size() != edf.timeline.num_total_epochs() ) 
    Helper::halt( "epoch annotation file " + f + " contains " 
		  + Helper::int2str( (int)a.size() ) + " epochs but expecting " 
		  + Helper::int2str( edf.timeline.num_total_epochs() ) );


  //
  // create and store proper annotation events
  //
  
  annot_t::map_epoch_annotations( edf , 
				  a , 
				  f , 
				  edf.timeline.epoch_len_tp() , 
				  edf.timeline.epoch_increment_tp() );
  
}


// ANNOT-MASK : add (internally) a MASK corresponding to included (or excluded epochs)
void proc_annot_mask( edf_t & edf , param_t & param )
{
  // default annot name = "E"
  const std::string tag = param.has( "inc" ) ? param.value( "inc" ) : "E" ;
  edf.timeline.add_mask_annot( tag );
}

// DUMP-MASK : output the current mask as an .annot file
void proc_dump_mask( edf_t & edf , param_t & param )
{
  
  edf.timeline.dumpmask( param );
  
  // nb. dumpmask() now allows adding annots +/- output
  // // otherwise, create an ANNOT file from the MASK, i.e. for viewing
  // std::string tag = param.requires( "tag" );
  // std::string path = param.has( "path" ) ? param.value("path") : ".";
  // bool no_id = ! param.has( "no-id" );
  // edf.timeline.mask2annot( path, tag , no_id ); 
}

// CHEP : dump, or convert from CHEP->MASK
// in edf/chep.cpp


// COUNT-ANNOTS : show all annotations for the EDF // REDUNDANT
void proc_list_annots( edf_t & edf , param_t & param )
{
  summarize_annotations( edf , param );
}


// ESPAN : per epoch, give # secs/proportion spanned by annotations
void proc_espan(  edf_t & edf , param_t & param )
{
  edf.annotations->espan( edf, param );
}

// MAKE-ANNOTS : make new annotations based on pairwse intersection,union, overlap, etc
void proc_make_annots( edf_t & edf , param_t & param )
{
  edf.annotations->make( param , edf );
}


// META : set annotation meta-data
void proc_set_annot_metadata( edf_t & edf , param_t & param )
{
  edf.timeline.set_annot_metadata( param );
}

// WRITE-ANNOTS : write all annots to disk
void proc_write_annots( edf_t & edf , param_t & param )
{
  if ( param.has( "annot-dir" ) )
    {
      if ( param.has( "file" ) ) Helper::halt( "cannot specify both annot-dir and file" );
      edf.annotations->write( "__special_make_dir__" , param , edf );
    }
  else
    edf.annotations->write( param.requires( "file" ) , param , edf );
}

// EXTEND : make single point annots longer
void proc_extend_annots( edf_t & edf , param_t & param )
{
  edf.annotations->extend( param );
}


// ANNOTATE : annotate one annot based on other annotation(s)
//  (for this one EDF) - also see command-line variant to read
//  in annots from multiple individuals for a joint test
void proc_annotate( edf_t & edf , param_t & param )
{
  annotate_t annotate( edf , param );
}


// AXA : annotation cross-tabs
void proc_annot_crosstabs( edf_t & edf , param_t & param )
{
  edf.timeline.annot_crosstabs( param );
}

// A2S : make signbal from ANNOTS
void proc_annot2signal( edf_t & edf , param_t & param )
{
  edf.timeline.annot2signal( param );
}

// S2A : make annot from a signal
void proc_signal2annot( edf_t & edf , param_t & param )
{
  edf.timeline.signal2annot( param );
}


// A2C : make a cache from an annotation
void proc_annot2cache( edf_t & edf , param_t & param )
{
  edf.timeline.annot2cache( param );
}

// C2A : make an annotation from a cache
void proc_cache2annot( edf_t & edf , param_t & param )
{
  edf.timeline.cache2annot( param );
}

// MEANS : signal means conditional on annotations
void proc_sig_annot_mean( edf_t & edf , param_t & param )
{
  edf.timeline.signal_means_by_annot( param );
}

// TABULATE : assume discrete values for a signal, and get counts
void proc_sig_tabulate( edf_t & edf , param_t & param )
{
  edf.tabulate( param );
}

// ANNOTS : list all annotations
void proc_list_all_annots( edf_t & edf , param_t & param )
{
  edf.timeline.list_all_annotations( param );
}


// ANNOTS-SPANNING : list all annotations
void proc_list_spanning_annots( edf_t & edf , param_t & param )
{
  edf.timeline.list_spanning_annotations( param );
}

// TIME-TRACK : make EDF+
void proc_timetrack( edf_t & edf , param_t & param )
{
  edf.add_time_track();
}

// RESTRUCTURE : flush masked records
void proc_restructure( edf_t & edf , param_t & param )
{
  // just drop MASK'ed records, then reset mask
  const bool FORCE_RESTRUCTURE = false;
  const bool VERBOSE_OUTPUT = param.has( "verbose" );
  const bool PRESERVE_CACHE = param.has( "preserve-cache" ) ? param.yesno( "preserve-cache" ) : false ; 
  edf.restructure( FORCE_RESTRUCTURE , VERBOSE_OUTPUT , PRESERVE_CACHE );

  // optionally, require a certain number of epochs to be left?
  if ( param.has( "require" ) )
    {
      const int r = param.requires_int( "require" );
      if ( r > 0 ) 
	{
	  edf.timeline.ensure_epoched();
	  
	  if ( edf.timeline.num_epochs() < r ) 
	    {         
	      logger << " ** warning: after RESTRUCTURE: "
		     << "required " << r << " epochs "
		     << "but observed only " << edf.timeline.num_epochs()  << " (setting empty flag)\n";
	      globals::empty = true;
	    }
	  return;
	}
    }
}


// DUMP-RECORDS : show all records (i.e. raw data)
void proc_record_dump( edf_t & edf , param_t & param )
{
  edf.add_time_track();  
  edf.record_dumper( param );
}

// SEGMENTS : show all contiguous segments
// (and optionally, add annotations to this effect)
void proc_dump_segs( edf_t & edf , param_t & param )
{
  edf.seg_dumper( param );
}


// RECS : simple table of records, epochs
void proc_record_table( edf_t & edf , param_t & param )
{
  edf.record_table( param );
}

// STAGE : set and display sleep stage labels (verbose = F)
// STAGE : + eannot=<file> option --> write as .eannot
// HYPNO : verbose report on sleep STAGES     (verbose = T)
void proc_sleep_stage( edf_t & edf , param_t & param , bool verbose )
{
  
  std::string wake   = param.has( "W" )  ? param.value("W")  : "" ; 
  std::string nrem1  = param.has( "N1" ) ? param.value("N1") : "" ; 
  std::string nrem2  = param.has( "N2" ) ? param.value("N2") : "" ; 
  std::string nrem3  = param.has( "N3" ) ? param.value("N3") : "" ; 
  std::string nrem4  = param.has( "N4" ) ? param.value("N4") : "" ; 
  std::string rem    = param.has( "R" )  ? param.value("R")  : "" ;
  std::string lights = param.has( "L" )  ? param.value("L")  : "" ; 
  std::string misc   = param.has( "?" )  ? param.value("?")  : "" ; 
  bool force_remake  = param.has( "force" );
  
  std::string eannot = param.has( "eannot" ) ? param.value( "eannot" ) : "" ;
  if ( eannot != "" && verbose ) Helper::halt( "cannot use eannot with HYPNO" );

  // simmple dump to standard out for STAGE
  if ( param.has( "min" ) ) eannot = "."; // code --> std::cout
  
  // either read these from a file, or display
  
  if ( param.has( "file" ) )
    {
      std::vector<std::string> ss = Helper::file2strvector( param.value( "file" ) );
      bool okay = edf.timeline.hypnogram.construct( &edf.timeline , param , verbose , ss );
      if ( ! okay ) return; // i.e. if no valid annotations found
    }
  else
    {      
      bool okay = edf.annotations->make_sleep_stage( edf.timeline, force_remake,
							     wake , nrem1 , nrem2 , nrem3 ,
							     nrem4 , rem , lights, misc );
      if ( ! okay ) return; // e.g. overlapping stages
      okay = edf.timeline.hypnogram.construct( &edf.timeline , param , verbose ); 
      if ( ! okay ) return; // i.e. if no valid annotations found
    }

  // epoch level output for HYPNO?
  bool epoch_lvl_output =param.has( "epoch" );

  // cycles and transitions?
  bool cycle_lvl_output = param.has( "verbose" ) ? param.yesno( "verbose" ) : true ;
    
  // optionally, add annotations
  // allow old annot-cycles form
  
  std::string annot_prefix = "";

  if ( param.has( "annot" ) )
    {
      if ( param.empty( "annot" ) )
	annot_prefix = "h";
      else
	annot_prefix = param.value( "annot" );
    }
  else if ( param.has( "annot-cycles" ) )
    {
      if ( param.empty( "annot-cycles" ) )
	annot_prefix = "h";
      else
	annot_prefix = param.value( "annot-cycles" );
    }
  
  
  // and output...

  edf.timeline.hypnogram.output( verbose , epoch_lvl_output , cycle_lvl_output, eannot , annot_prefix );

}


// ED : compute 'electrical distance' measure of bridging
void proc_elec_distance( edf_t & edf , param_t & param )
{
  dsptools::elec_distance( edf , param );
}

// L1OUT : leave-one-out validation via interpolation of all signals
void proc_leave_one_out( edf_t & edf , param_t & param )
{
  dsptools::leave_one_out( edf , param );
}

// INTERPOLATE : chep-based interpolation
void proc_chep_based_interpolation( edf_t & edf , param_t & param )
{
  dsptools::chep_based_interpolation( edf , param );
}

// SL : surface laplacian
void proc_surface_laplacian( edf_t & edf , param_t & param )
{
  dsptools::surface_laplacian_wrapper( edf , param );
}

// CLOCS : attach channel locations
void proc_attach_clocs( edf_t & edf , param_t & param )
{

  if ( ! param.has( "file" ) ) 
    {
      edf.clocs.set_default();
      return;
    }

  std::string filename = Helper::expand( param.requires( "file" ) );
  if ( ! Helper::fileExists( filename ) ) Helper::halt( "could not find " + filename );

  // assume a cartesian format
  edf.clocs.load_cart( filename , param.has( "verbose" ) );
}


// EMD : Empirical Mode Decomposition 
void proc_emd( edf_t & edf , param_t & param )
{
  dsptools::emd_wrapper( edf , param );
}

// DFA : detrended fluctuation analysis
void proc_dfa( edf_t & edf , param_t & param )
{
  dsptools::dfa_wrapper( edf , param );
}

// ICA : fastICA on sample by channel matrix (whole trace)
void proc_ica( edf_t & edf , param_t & param )
{
  dsptools::ica_wrapper( edf , param );
}

// SVD : singular value decomposition on signals
void proc_svd( edf_t & edf , param_t & param )
{
  dsptools::svd_wrapper( edf , param );
}

// COH : calculate cross spectral coherence, using new/default code
void proc_coh( edf_t & edf , param_t & param )
{
  dsptools::coherence( edf , param );
}

// CORREL : correlation
void proc_correl( edf_t & edf , param_t & param )
{
  dsptools::correlate_channels( edf , param );
}

// ACF : autocorrelation function
void proc_acf( edf_t & edf , param_t & param )
{
  dsptools::autocorr_channels( edf , param );
}

// MI : mutual information
void proc_mi( edf_t & edf , param_t & param )
{
  dsptools::compute_mi( edf , param );
}

// CC : gerneral connectivity and coupling metrics, using wavelets or filter-Hilbert
void proc_conncoupl( edf_t & edf , param_t & param )
{
  dsptools::connectivity_coupling( edf , param );
}

// SHIFT : shift one or more signals by X samples
void proc_shift( edf_t & edf , param_t & param )
{
  dsptools::shift( edf , param );
}

// SCRAMBLE : randomly scramble a signal (complete perm)
void proc_scramble( edf_t & edf , param_t & param )
{
  dsptools::scramble( edf , param );
}


// CACHE : internal command to dump cache contents (debugging)
void proc_dump_cache( edf_t & edf , param_t & param )
{

  if ( param.has( "test" ) )
    {
      ctest2( edf );
      return;
    }

  // set up to record from output stream
  if ( param.has( "record" ) )
    {
      // record=command,variable,{strata}
      std::vector<std::string> tok = param.strvector( "record" );
      if ( tok.size() < 2 ) Helper::halt( "record=command,variable,{strata}" );

      std::set<std::string> fac;
      for (int i=2; i<tok.size(); i++) fac.insert( tok[i] );

      // assume numeric
      bool is_text = param.has( "text" ) || param.has( "str" );
      //      bool is_int  = param.has( "int" ) || param.has( "integer" );

      const std::string cname = param.requires( "cache" );
      const std::string cmd_name = Helper::toupper( tok[0] );
      const std::string var_name = Helper::toupper( tok[1] );
      const std::string fac_name = Helper::toupper( Helper::stringize( fac ) );

      if ( is_text )
	{
	  // makes if does not exist
	  cache_t<std::string> * cache = edf.timeline.cache.find_str( cname );
	  writer.cache( cache );
	  writer.cache_string( cmd_name , var_name, fac );
	  if ( param.has( "uppercase-keys" ) ) cache->set_uppercase_keys( true ); 
	}
      // else if ( is_int )
      // 	{
      // 	  cache_t<int> * cache = edf.timeline.cache.find_int( cname );
      // 	  writer.cache( cache );
      // 	  writer.cache_integer( cmd_name , var_name , fac );
      // 	}
      else
	{
	  cache_t<double> * cache = edf.timeline.cache.find_num( cname );
	  writer.cache( cache );
	  writer.cache_numeric( cmd_name , var_name , fac );
	  if ( param.has( "uppercase-keys" ) ) cache->set_uppercase_keys( true ); 
	}
      
      logger << "  caching output from " << cmd_name << ", variable = " << var_name ;
      if ( fac.size() != 0 ) logger << " (strata = " << fac_name << ")";
      logger << " to cache " << cname << "\n";
								       
    }
  
  // clear / load / dump 

  if ( param.has( "clear" ) )
    {
      // clear all caches
      edf.timeline.cache.clear();
    }

  // load from prior CACHE save 

  if ( param.has( "load" ) )
    {
      const std::string & filename = param.value( "load" );
      if ( ! Helper::fileExists( filename ) ) 
	Helper::halt( "cannot find " + filename );
      edf.timeline.cache.load( filename );
    }

  // load from long-format output (destrat out.db +COMMAND -r F1 F2)
  // expects ID, factors, variables... 
  // only load from that individual
  
  if ( param.has( "import" ) ) 
    {
      const std::string & filename = param.value( "import" );

      if ( ! Helper::fileExists( filename ) )
	Helper::halt( "cannot find " + filename );
      
      std::set<std::string> factors;

      if ( param.has( "factors" ) ) 
	factors = param.strset( "factors" );

      std::set<std::string> vars;
      if ( param.has( "v" ) ) 
	vars = param.strset( "v" );
      
      edf.timeline.cache.import( filename , 
				 param.requires( "cache" ) , 				 
				 edf.id , 
				 factors ,
				 param.has( "v" ) ? &vars : NULL );
      
    }

  
  if ( param.has( "dump" ) )
    {

      // cache types: int, num and tp
      int int_cache = param.has( "int" );
      int str_cache = param.has( "str" );
      int num_cache = param.has( "num" );
      int tp_cache = param.has( "tp" );
      
      if ( str_cache + int_cache + num_cache + tp_cache != 1 )
	Helper::halt( "need to specify one of int, str, num or tp cache types" );
      
      std::string cname;
      if ( int_cache ) cname = param.value( "int" );
      else if ( str_cache ) cname = param.value( "str" );
      else if ( num_cache ) cname = param.value( "num" );
      else cname = param.value( "tp" );
      
      if ( int_cache )
	{
	  cache_t<int> * cache = edf.timeline.cache.find_int( cname );
	  if ( cache == NULL ) Helper::halt( "could not find int-cache " + cname );
	  std::cout << "cache: " << cname << "[int]\n";
	  std::cout << cache->print();
	}
      else if ( str_cache )
	{
	  cache_t<std::string> * cache = edf.timeline.cache.find_str( cname );
	  if ( cache == NULL ) Helper::halt( "could not find str-cache " + cname );
	  std::cout << "cache: " << cname << "[str]\n";
	  std::cout << cache->print();
	}  
      else if ( num_cache )
	{
	  cache_t<double> * cache = edf.timeline.cache.find_num( cname );
	  if ( cache == NULL ) Helper::halt( "could not find num-cache " + cname );
	  std::cout << "cache: " << cname << "[num]\n";
	  std::cout << cache->print();	  
	}
      else if ( tp_cache )
	{
	  cache_t<uint64_t> * cache = edf.timeline.cache.find_tp( cname );
	  if ( cache == NULL ) Helper::halt( "could not find tp-cache " + cname );
	  std::cout << "cache: " << cname << "[tp]\n";
	  std::cout << cache->print();
	}
      
    }
}


// SIGGEN : add/generate artificial signals to 1+ channels
void proc_siggen( edf_t & edf , param_t & param )
{
  dsptools::siggen( edf, param );
}

// SIMUL : simulate a time series from a PSD
void proc_simul( edf_t & edf , param_t & param )
{
  dsptools::simul( edf , param );
}

// SPIKE : spike in a new bit of signal
void proc_spike( edf_t & edf , param_t & param )
{

  // create a new signal?
  std::string ns = "";
  
  if ( param.has( "new" ) ) ns = param.value( "new" );
  
  signal_list_t from_signal = edf.header.signal_list( param.requires( "from" ) );  
  signal_list_t to_signal   = edf.header.signal_list( param.requires( "to" ) );  
  
  if ( from_signal.size() != 1 ) Helper::halt( "no from={signal}" );
  if ( to_signal.size() != 1 ) Helper::halt( "no to={signal}" );
  
  int s1 = to_signal(0);
  int s2 = from_signal(0);
  
  double wgt = param.requires_dbl( "wgt" );

  spike_signal( edf , s1 , s2 , wgt , ns );
}

// PAC : phase amplitude coupling
void proc_pac( edf_t & edf , param_t & param )
{
  dsptools::pac( edf , param );
}


// GED : generalized eigendecomposition
void proc_ged( edf_t & edf , param_t & param )
{
  ged_wrapper( edf, param );
}

// CFC : generic cross-frequency coupling methods (other than PAC as above)
void proc_cfc( edf_t & edf , param_t & param )
{
  dsptools::cfc( edf , param );
}

// HB : Hypoxic burden
void proc_hypoxic_burden( edf_t & edf , param_t & param )
{
  hb_t hb( edf , param );
}

// SUPPRESS-ECG : ECG supression
void proc_ecgsuppression( edf_t & edf , param_t & param )
{
  dsptools::ecgsuppression( edf , param );
}

// HRV : heart rate variability metrics
void proc_hrv( edf_t & edf , param_t & param )
{
  dsptools::hrv( edf , param );
}

// BPM : get beats per min from ECG
void proc_bpm( edf_t & edf , param_t & param )
{
  dsptools::bpm( edf , param );
}

// ORDER : set order of signals ( in internal EDF)
void proc_order_signals( edf_t & edf , param_t & param )
{
  // retain only specified set (all must exist)
  // swap in EDF internal
  edf.set_order( param );
}


// READ: just read all data in, do nothing else (for timing)
void proc_read_signal( edf_t & edf , param_t & param )
{
  edf.preread(param);
}

// COPY : mirror a signal
void proc_copy_signal( edf_t & edf , param_t & param )
{
  
  signal_list_t originals = edf.header.signal_list( param.requires( "sig" ) );
  
  const bool newlab = param.has( "new" );
  const bool pretag = param.has( "pretag" );
  const bool posttag = param.has( "tag" );

  if ( (int)newlab + (int)pretag + (int)posttag != 1 )
    Helper::halt( "requires one of tag, pretag or new" );
  
  const std::string tag = posttag ? param.value( "tag" ) : ( pretag ? param.value( "pretag" ) : param.value( "new" ) );

  if ( newlab && originals.size() > 1 ) Helper::halt( "cannot specify 'new' with more than one 'sig'" );
  
  for (int s=0;s<originals.size();s++)
    {

      if ( edf.header.is_data_channel( originals(s) ) )
	{
	  const std::string new_label = newlab ? tag :
	    ( posttag ? originals.label( s ) + tag : tag + originals.label( s ) );
	  
	  if ( ! edf.header.has_signal( new_label ) )
	    {
	      logger << "  copying " << originals.label(s) << " to " << new_label << "\n";
	      edf.copy_signal( originals.label(s) , new_label );
	    }
	  else
	    logger << "  *** " << new_label << " already exists, skipping copying operation\n";
	}
    }
}

// ENFORCE-SR : drop/alter signals based on SR requirements (for record size)
void proc_enforce_signals( edf_t & edf , param_t & param )
{
  
  // to enable clean record-size conversion, this first drops any signals which would not
  // be represented by a N-second record size (i.e. requires integer Hz sample rate)

  std::set<std::string> drops;
  
  const bool no_annotations = true; 

  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) , no_annotations );

  const int ns = signals.size();
  
  // default: 1-second EDF record size
  const double new_record_duration = param.has( "dur" ) ? param.requires_dbl( "dur" ) : 0 ; 
    
  std::vector<double> range;
  if ( param.has( "sr" ) )
    {
      range = param.dblvector( "sr" );
      if ( range.size() != 2 ) Helper::halt( "expecting sr=lwr,upr" ) ;
      if ( range[0] > range[1] ) Helper::halt( "expecting sr=lwr,upr" ) ;
    }

    if ( new_record_duration > 0 )
    logger << "  retaining channels that can be represented in an EDF record of " << new_record_duration << " second\n";
  if ( range.size() == 2 )
    logger << "  retaining channels with SR between " << range[0] << " and " << range[1] << "\n";

  for (int s=0; s<ns; s++)
    {      
      int    nsamples = edf.header.n_samples[ signals(s) ];      
      double fs = (double)nsamples / edf.header.record_duration;

      // does new record size contain an integer number of sample points?
      if ( new_record_duration > 0 )
	{
	  double implied = new_record_duration * fs;
	  int   new_nsamples1 = implied;	  
	  if ( fabs( (double)new_nsamples1 - implied ) > 0 )
	    drops.insert( signals.label(s) );
	}
      
      // drop based on whether SR is within range?
      if ( range.size() == 2 )
	{
	  if ( fs < range[0] || fs > range[1] )
	    drops.insert( signals.label(s) );
	}
    }

  //
  // drop channels as needed
  //
  
  if ( drops.size() > 0 ) logger << "  dropping channels:";
  std::set<std::string>::const_iterator dd = drops.begin();
  while ( dd != drops.end() )
    {
      if ( edf.header.has_signal( *dd ) )
	{	  	  
	  logger << " " << *dd ;
	  int s = edf.header.signal( *dd );
	  edf.drop_signal( s );	  
	}
	++dd;
    }
  if ( drops.size() > 0 ) logger << "\n";
 
}


// RENAME : rename signals
void proc_rename( edf_t & edf , param_t & param )
{

  // either from a file, or the command line
  if ( param.has( "file" ) )
    {
      if ( param.has( "new" ) ) 
	Helper::halt( "cannot specify both file and sig/new" );
      
      std::vector<std::string> old_signals, new_signals;
      std::set<std::string> newset;

      const std::string fname = Helper::expand( param.value( "file" ) );
      if ( ! Helper::fileExists( fname ) )
	Helper::halt( "could not open " + fname );
      
      std::ifstream I1( fname.c_str() , std::ios::in );
      while ( ! I1.eof() )
	{
	  std::string line;
	  Helper::safe_getline( I1 , line);
	  if ( I1.eof() || line == "" ) continue;	  
	  std::vector<std::string> tok2 = Helper::parse( line , "\t" );
	  if ( tok2.size() != 2 ) 
	    Helper::halt( "expecting two tab-delimited values: " + line );

	  const std::string s1 = tok2[0];
	  const std::string s2 = tok2[1];
	  	  
	  const bool old_exists = edf.header.has_signal( s1 );
	  const bool new_exists = edf.header.has_signal( s2 );
	  
	  // new mappings must be unique
	  if ( new_exists )
	    Helper::halt( "'new' signal labels cannot already exist in the EDF" );
	  
	  // just ignore if original channel does not exist
	  if ( old_exists )
	    {
	      old_signals.push_back( s1 );
	      new_signals.push_back( s2 );
	      newset.insert( s2 );
	    }
	  
	  // next pair
	}

      if ( newset.size() != new_signals.size() )
	Helper::halt( "cannot have duplicate labels in new" );
      
      for (int s=0; s<old_signals.size(); s++)
	{
	  logger << "  renaming [" << old_signals[s] << "] as [" << new_signals[s]  << "]\n";
	  edf.header.rename_channel( old_signals[s] , new_signals[s] );
	}
      
      return;
    }
  
  //
  // Otherwise, take input from command line
  //

  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );
  std::vector<std::string> new_signals = param.strvector( "new" );
  
  if ( signals.size() != new_signals.size() )
    Helper::halt( "number of channels for 'sig' and 'new' must match" );

  const int ns = signals.size();
  
  // check that all new labels are in fact new
  std::set<std::string> newset;
  for (int s=0; s<ns; s++)
    {
      if ( edf.header.has_signal( new_signals[s] ) )
	Helper::halt( "'new' signal labels cannot already exist in the EDF" );
      newset.insert( new_signals[s] );
    }
  
  if ( newset.size() != new_signals.size() )
    Helper::halt( "cannot have duplicate labels in new" );
  
  for (int s=0; s<ns; s++)
    {
      logger << "  renaming [" << signals.label(s) << "] as [" << new_signals[s]  << "]\n";
      edf.header.rename_channel(  signals.label(s) , new_signals[s] );
    }

}

// DROP-ANNOTS : drop one or more annotations
void proc_drop_annots( edf_t & edf , param_t & param )
{
  
  std::vector<std::string> drops;

  if      ( param.has( "annot" ) ) drops = param.strvector( "annot" ) ;
  else if ( param.has( "annots" ) ) drops = param.strvector( "annots" ) ;
  else edf.annotations->drop();
      
  if ( drops.size() > 0 ) 
    edf.annotations->drop( &drops );
  
}
  

// SIGNALS : drop one or more signal
void proc_drop_signals( edf_t & edf , param_t & param )
{
  
  std::set<std::string> keeps, drops;
  std::vector<std::string> picks; // order matters here

  if ( param.has( "keep" ) ) keeps = param.strset( "keep" );

  if ( param.has( "keep" ) && param.has( "req" ) ) 
    Helper::halt( "cannot specify both keep and req" );

  bool req = param.has( "req" ) ;
  if ( param.has( "req" ) ) keeps = param.strset( "req" );

  bool pick = param.has( "pick" );
  if ( pick && req ) Helper::halt( "cannot specify pick and req together" );
  if ( pick && param.has("drop") ) Helper::halt( "cannot specify pick and drop together" );
  if ( pick && param.has("keep") ) Helper::halt( "cannot specify pick and keep together" );
  if ( pick ) picks = param.strvector( "pick" );
  std::string pick_choice = "";

  std::string pick_rename = param.has( "rename" ) ? param.value( "rename" ) : "" ; 
  if ( edf.header.has_signal( pick_rename ) )
    Helper::halt( "rename choice already exists" );
  
  if ( param.has( "drop" ) ) drops = param.strset( "drop" );
  
  if ( param.has( "keep" ) && param.has( "drop" ) )
    Helper::halt( "can only specify keep or drop with SIGNALS" );
  
  if ( ! ( param.has( "pick" ) || param.has( "keep" ) || param.has( "drop" ) || param.has( "req" ) ) ) 
    Helper::halt( "need to specify keep, drop, pick or req with SIGNALS" );

  
  
  //
  // pick list?  use this to define a drop list
  //

  if ( picks.size() > 0 )
    {
      bool picked = false; 
      for (int p=0; p<picks.size(); p++)
	{
	  if ( edf.header.has_signal( picks[p] ) )
	    {
	      if ( ! picked )
		{
		  logger << "  picked " << picks[p] << "\n";
		  picked = true;
		  pick_choice = picks[p];
		}
	      else // add to the drop list
		{
		  drops.insert( picks[p] );
		}
	    }
	}
    }
  
  // if a keep list is specified, means we keep 
  if ( keeps.size() > 0 )
    {

      //
      // check that any required signals are present
      //

      if ( req ) 
	{
	  std::set<std::string>::const_iterator ss = keeps.begin();
	  while ( ss != keeps.end() )
	    {
	      if ( ! edf.header.has_signal( *ss ) )
		{
		  logger << "  *** could not find requested signal: " << *ss << "\n";
		  logger << "  *** quitting for this individual\n";
		  globals::problem = true;
		  return; 
		}
	      ++ss;
	    }
	}

      //
      // Otherwise, keep what we can 
      //

      const int ns = edf.header.ns;

      for (int s = 0 ; s < ns ; s++ )
	{
	  std::string label = edf.header.label[s];

	  // is this signal on the keep list?
	  if ( keeps.find( label ) == keeps.end() )
	    {
	      
	      // or, does this signal have an alias that is on the keep list?
	      if ( cmd_t::label_aliases.find( label ) != cmd_t::label_aliases.end() )
		{
		  //std::cout << " has alias " << cmd_t::label_aliases[ label ]  << "\n";
		  if ( keeps.find( cmd_t::label_aliases[ label ] ) == keeps.end() )
		    {
		      //std::cout << "drps " << label << "\n";
		      drops.insert( label );
		      // OR ?? drops.insert( cmd_t::label_aliases[ label ] ) ; should be equiv.
		    }
		}
	      else
		drops.insert( label );
	    }
	  
	} // next signal
    }


  if ( drops.size() > 0 ) logger << "  dropping channels:";
  std::set<std::string>::const_iterator dd = drops.begin();
  while ( dd != drops.end() )
    {
      if ( edf.header.has_signal( *dd ) )
	{	  	  
	  logger << " " << *dd ;
	  int s = edf.header.signal( *dd );
	  edf.drop_signal( s );	  
	}
	++dd;
    }
  if ( drops.size() > 0 ) logger << "\n";

  //
  // rename picked channel to something else?
  //

  if ( pick_choice != "" && pick_rename != "" )
    {
      logger << "  renaming pick, from " << pick_choice << " to " << pick_rename << "\n";
      edf.header.rename_channel( pick_choice , pick_rename ) ; 
    }
  
}


// SLICE : pull out slices, based on 'file'
void proc_slice( edf_t & edf , param_t & param , int extract )
{
  
  // x = 1 (extract)
  // x = 0 (exclude)
  
  std::string filename = Helper::expand( param.requires( "file" ) );
  
  std::set<interval_t> intervals;

  if ( ! Helper::fileExists( filename ) ) Helper::halt( "could not find " + filename );


  std::ifstream IN1( filename.c_str() , std::ios::in );  
  while ( ! IN1.eof() )
    {
      interval_t interval;
      IN1 >> interval.start >> interval.stop;
      if ( IN1.eof() ) break;
      if ( interval.stop <= interval.start ) Helper::halt( "problem with interval line" );
      intervals.insert(interval);
    }
  IN1.close();

  logger << " read " << intervals.size() << " from " << filename << "\n";
  
  edf.slicer( intervals , param , extract );
  
}


// REMAP
void proc_remap_annots( edf_t & edf , param_t & param )
{
  // as if having originally 'remap' command, but apply these
  // after the fact, i.e. to already loaded/created annots

  if ( ! param.has( "file" ) ) Helper::halt( "requires file argument" );
  
  const std::vector<std::string> files = param.strvector( "file" );
  
  int remap_field = 0;
  if ( param.has( "remap-col" ) ) remap_field = 1;
  else if ( param.has( "optional-remap-col" ) ) remap_field = 2;

  // be default, allow spaces (i.e. for moonlight)
  const bool remap_spaces = param.has( "allow-spaces" ) ? param.yesno( "allow-spaces" ) : false;
  const bool remap_verbose = param.has( "verbose" );
  
  int mapped = edf.annotations->remap( files , remap_field , remap_spaces , remap_verbose );

  logger << "  remapped " << mapped << " annotations\n";
  
}


// CANONICAL
void proc_canonical( edf_t & edf , param_t & param )
{
  
  if ( ! param.has( "legacy" ) )
    {
      canonical_t canonical( edf , param );  
      return;
    }
  
  // dry-run or make the actual signals?
  bool make_signals = ! param.has( "check" );

  // just try to guess... no file specification
  if ( param.has( "guess" ) )
    {
      edf.guess_canonicals( param , make_signals );
      return;
    }

  // canonical signal file
  if ( ! ( param.has( "file" ) || param.has( "files" ) ) )
    Helper::halt( "one or more definition files required, file=cs1.txt,cs2.txt" );
  
  const std::vector<std::string> files = param.strvector( param.has( "file" ) ? "file" : "files" );
  
  // (optional) group for the canonical file?
  const std::string group = param.has( "group" ) ? param.value( "group" ) : "." ; 

  // add prefix to canonical labels? 
  const std::string prefix = param.has( "prefix" ) ? param.value( "prefix" ) : "" ;

  // drop all non-canonical signals from EDF after processing?
  const bool drop_originals = param.has( "drop-originals" );

  // if ( drop_originals && ! make_signals )
  //   Helper::halt( "cannot have drop-originals and check options together for CANONICAL" );

  // cs = additional subset of canonical signals to focus on
  if ( ! param.has( "cs" ) )    
    cansigs_t cs0 = edf.make_canonicals( files , group , make_signals , drop_originals , prefix );
  else
    {
      const std::set<std::string> cs = param.strset( "cs" );
      cansigs_t cs0 = edf.make_canonicals( files , group , make_signals , drop_originals , prefix , &cs );
    }
}


// ADJUST: Adjust signals by ICs 
void proc_adjust( edf_t & edf , param_t & param )
{
  dsptools::ica_adjust( edf , param );
}

// REFERENCE : re-reference tracks
void proc_reference( edf_t & edf , param_t & param )
{
  std::string sigstr = param.requires( "sig" );
  signal_list_t signals = edf.header.signal_list( sigstr );

  signal_list_t references;
  std::string refstr = param.requires( "ref" );
  if ( refstr != "." ) references = edf.header.signal_list( refstr );

  // if new channel label given, then also rename signal 
  // REFERENCE sig=C3 ref=A2 new=C3_A2
  //   this leaves C3 as is, and make a new channel called 'C3_A2'

  bool make_new = param.has( "new" );

  std::vector<std::string> new_channels;
  
  const bool pairwise = param.has( "pairwise" ) ;
  if ( make_new )
    {
      new_channels = param.strvector( "new" );
      if ( ! pairwise )
        if ( new_channels.size() != 1 )
          Helper::halt( "expecting a single label for new" );
    }
  
  int new_sr = 0;
  if ( make_new && param.has( "sr" ) ) new_sr = param.requires_int( "sr" );
  
  if ( pairwise )
    edf.pairwise_reference( signals , references , make_new , new_channels , new_sr , false );
  else
    edf.reference( signals , references , make_new , new_channels[0] , new_sr , false );
  
}

// Remove reference
void proc_dereference( edf_t & edf , param_t & param )
{
  std::string sigstr = param.requires( "sig" );
  signal_list_t signals = edf.header.signal_list( sigstr );
  
  signal_list_t references;
  std::string refstr = param.requires( "ref" );
  if ( refstr != "." ) references = edf.header.signal_list( refstr );

  bool make_new = param.has( "new" );
  std::vector<std::string> new_channels;
  
  const bool pairwise = param.has( "pairwise" ) ;
  if ( make_new )
    {
      new_channels = param.strvector( "new" );
      if ( ! pairwise )
	if ( new_channels.size() != 1 )
	  Helper::halt( "expecting a single label for new" );
    }
  
  int new_sr = 0;
  if ( make_new && param.has( "sr" ) )
    new_sr = param.requires_int( "sr" );

  if ( pairwise ) 
    edf.pairwise_reference( signals , references , make_new , new_channels , new_sr , true );
  else
    edf.reference( signals , references , make_new , new_channels[0] , new_sr , true );
  
}


// RECSIZE: change record size for one or more signals
void proc_rerecord( edf_t & edf , param_t & param )
{
  double rs = param.requires_dbl( "dur" ); 
  
  logger << "  altering record size from " << edf.header.record_duration << " to " <<  rs << " seconds\n";
  
  edf.reset_record_size( rs );
  
  logger << "  *** now WRITE'ing EDF to disk, and will set 'problem' flag to skip to next EDF *** \n";
  
  proc_write( edf , param );
  globals::problem = true;
}


// SCALE
void proc_setscale( edf_t & edf , param_t & param )
{
  
  const bool NO_ANNOTS = true;
  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) , NO_ANNOTS );

  // expects min,max
  std::vector<double> minmax ;
  const bool do_minmax = param.has( "min-max" ) ;
  if ( do_minmax )
    {
      minmax = param.dblvector( "min-max" );
      if ( minmax.size() != 2 || minmax[0] >= minmax[1] )
	Helper::halt( "expecting two valies max-max=a,b  where a is lower than b" );
    }
  
  // clip also?
  double clip_min , clip_max;
  bool   do_clip_min = false , do_clip_max = false;

  if ( param.has( "clip-min" ) ) 
    {
      do_clip_min = true;
      clip_min = param.requires_dbl( "clip-min" );      
    }

  if ( param.has( "clip-max" ) ) 
    {
      do_clip_max = true;
      clip_max = param.requires_dbl( "clip-max" );      
    }
  
  if ( ! ( do_minmax | do_clip_min | do_clip_max ) ) 
    {
      logger << "  nothing to do, returning\n";
      return;
    }
  
  const int ns = signals.size();
  for (int s=0;s<ns;s++)
    {
      if ( minmax.size() == 2 ) 
	edf.set_scale( signals(s) , &(minmax[0]) , &(minmax[1]) , do_clip_min ? &clip_min : NULL , do_clip_max ? &clip_max : NULL );
      else
	edf.set_scale( signals(s) , NULL , NULL, do_clip_min ? &clip_min : NULL , do_clip_max ? &clip_max : NULL );
    }  
  
}


// uV or mV : set units for tracks
void proc_scale( edf_t & edf , param_t & param , const std::string & sc )
{
  std::string sigstr = param.requires( "sig" );
  signal_list_t signals = edf.header.signal_list( sigstr );
  const int ns = signals.size();  
  for (int s=0;s<ns;s++) 
    edf.rescale( signals(s) , sc );
}

// MINMAX : set EDF header to have identical min/max (physical) values
void proc_minmax( edf_t & edf , param_t & param )
{
  std::string sigstr = param.requires( "sig" );
  signal_list_t signals = edf.header.signal_list( sigstr );

  // clip at these values? (and set PMIN/PMAX for these signals)

  const bool clip_min = param.has( "min" ) ;
  const bool clip_max = param.has( "max" );

  // by default, if |value| smaller than specified, leave as is
  // (i.e. only clip, do not expand range) ; if ! 
  
  const bool force = param.has( "force" );

  // if not 'clip' then MINMAX sets all PMIN/PMAX to be the same
  if ( ! ( clip_min || clip_max ) ) 
    edf.minmax( signals );
  else 
    {
      const double pmin = clip_min ? param.requires_dbl( "min" ) : 0 ;
      const double pmax = clip_max ? param.requires_dbl( "max" ) : 0 ;      
      edf.minmax( signals , clip_min ? &pmin : NULL , clip_max ? &pmax : NULL , force );      
    }
}

// ROLLING-NORM
void proc_rolling_norm( edf_t & edf , param_t & param )
{
  dsptools::rolling_standardize( edf, param );
}

// STANDARDIZE : robust winsorization and norming for each signal (done per whole signal or epoch) 
void proc_standardize( edf_t & edf , param_t & param )
{
  dsptools::standardize( edf , param );
}

// RAI : REM atonia index
void proc_rai( edf_t & edf , param_t & param  )
{
  dsptools::rai( edf, param );
}

// CLIP
void proc_clip( edf_t & edf, param_t & param )
{
  dsptools::clip( edf, param );
}

// RECTIFY 
void proc_rectify( edf_t & edf , param_t & param  )
{
  dsptools::rectify( edf , param );
}

// FLIP : change polarity of signal
void proc_flip( edf_t & edf , param_t & param  )
{
  std::string sigstr = param.requires( "sig" );
  signal_list_t signals = edf.header.signal_list( sigstr );
  const int ns = signals.size();  

  // track which signals are flipped
  for (int s=0;s<ns;s++) 
    {
      writer.level( signals.label(s) , globals::signal_strat );
      writer.value( "FLIP" , 1 );
      edf.flip( signals(s) );
    }
  writer.unlevel( globals::signal_strat );
}

// REVERSE : reverse signal in time domain
void proc_reverse( edf_t & edf , param_t & param  )
{
  std::string sigstr = param.requires( "sig" );
  signal_list_t signals = edf.header.signal_list( sigstr );
  const int ns = signals.size();  

  // track which signals are flipped
  for (int s=0;s<ns;s++) 
    {
      writer.level( signals.label(s) , globals::signal_strat );
      writer.value( "REVERSE" , 1 );
      edf.reverse( signals(s) );
    }
  writer.unlevel( globals::signal_strat );
}



//
// Parse special command-line arguments (e.g. path, silent, etc) 
//
	      
void cmd_t::parse_special( const std::string & tok0 , const std::string & tok1 )
{

  // always record as a standard variable also (i.e. for lookup in
  // lunapi, etc) n.b. previously, we'd skip this and only add
  // non-special vars

  // ... except handle the special case of 'sig' which _appends_ to an
  // existing list unless '.'

  
  


  
  //
  // ------------------------------------------------------------
  //

  if ( tok0 == "sig" )
    {
      std::string curr = cmd_t::vars[ tok0 ];
      // wipe?
      if ( tok1 == "." || tok1 == "" )
	cmd_t::vars.erase( cmd_t::vars.find( "sig" ) );
      else // add to list
	{
	  if ( curr == "." || curr == "" )
	    cmd_t::vars[ tok0 ] = tok1 ;
	  else
	    cmd_t::vars[ tok0 ] = curr + "," + tok1 ;
	}
    }
  else // just add as usual
    cmd_t::vars[ tok0 ] = tok1;


  //
  // now process sig opt via cmd_t::signallist
  //
  
  if ( Helper::iequals( tok0 , "sig" ) )
    {
      // wipe?
      if ( tok1 == "." || tok1 == "" )
	cmd_t::signallist.clear();
      else // or append to list?
	{
	  std::vector<std::string> tok2 = Helper::quoted_parse( tok1 , "," );		        
	  for (int s=0;s<tok2.size();s++) 
	    cmd_t::signallist.insert( globals::sanitize_everything ?
				      Helper::sanitize( Helper::unquote(tok2[s] ) ) :				  
				      Helper::unquote(tok2[s]) );
	}
      return;
    }
  

  //
  // now handle processing of all other special variables
  //

  // return signals alphabetically
  if ( Helper::iequals( tok0 , "order-signals" ) )
    {
      globals::order_signal_list_alphabetically = Helper::yesno( tok1 );
      return;
    }
  
  
  // no console logging?
  if ( Helper::iequals( tok0 , "silent" ) ) 
    {
      globals::silent = Helper::yesno( tok1 );
      return;      
    }

    
  // save console logging
  if ( Helper::iequals( tok0 , "log" ) )
    {
      globals::write_log = true;
      globals::log_file = Helper::expand( tok1 );
      return;
    }

      
  // allow to turn off log saving (e.g. via API)
  if ( Helper::iequals( tok0 , "write-log" ) )
    {
      globals::write_log = Helper::yesno( tok1 );      
      return;
    }

      
  if ( Helper::iequals( tok0 , "verbose" ) )
    {
      globals::verbose = Helper::yesno( tok1 );
      return;
    }

  

  if ( Helper::iequals( tok0 , "srand" ) )
    {
      
      uint64_t rng_seed = 0;
      //long unsigned 
      if ( ! Helper::str2int64( tok1  , &rng_seed ) )
	Helper::halt( "cannot set srand to " + tok1 );
      logger << "  setting RNG seed to " << rng_seed << "\n";
      CRandom::srand( rng_seed );
      return;
    }

  
  if ( Helper::iequals( tok0 , "devel" ) )
    {
      globals::devel = Helper::yesno( tok1 );
      return;
    }

  
  
  if ( Helper::iequals( tok0 , "legacy-hjorth" ) )
    {
      globals::legacy_hjorth = Helper::yesno( tok1 );
      return;
    }

  
  
  if ( Helper::iequals( tok0 , "show-assignments" ) )
    {
      globals::verbose_var_assignment = Helper::yesno( tok1 );
      return;
    }

  

  // mirror inputs in console log
  if ( Helper::iequals( tok0 , "mirror" ) )
    {
      globals::mirror = Helper::yesno( tok1 );
      return;
    }

  
  

  // debug mode: read digital values from EDF (i.e. do not translate
  // to physical values)

  if ( Helper::iequals( tok0 , "digital" ) )
    {
      globals::read_digital_values = Helper::yesno( tok1 );
      return;
    }

    
  // specify indiv(s) (i.e. can be used if ID is numeric)
  if ( Helper::iequals( tok0 , "id" ) )
    {
      globals::sample_list_ids = Helper::vec2set( Helper::parse( tok1 , "," ) );
      return;
    }

 

  // specify indiv(s) to skip (i.e. can be used if ID is numeric)
  if ( Helper::iequals( tok0 , "skip" ) )
    {
      globals::sample_list_ids_skips = Helper::vec2set( Helper::parse( tok1 , "," ) );
      return;
    }

  
  
  // do not read IDs in
  if ( Helper::iequals( tok0 , "anon" ) )
    {
      globals::anon = Helper::yesno( tok1 );
      return;
    }

  

  // read in US style annot dates (mm-dd-yy)
  if ( Helper::iequals( tok0 , "date-format" ) )
    {
      const std::string v = Helper::toupper( tok1 );
      if ( v == "MDY" ) globals::read_annot_date_format = MDY;
      else if ( v == "DMY" ) globals::read_annot_date_format = DMY;
      else if ( v == "YMD" ) globals::read_annot_date_format = YMD;
      return;
    } 

  
  // write in different date format
  if ( Helper::iequals( tok0 , "write-date-format" ) )
    {
      const std::string v = Helper::toupper( tok1 );
      if ( v == "MDY" )	globals::write_annot_date_format = MDY;
      else if ( v == "DMY" ) globals::write_annot_date_format = DMY;
      else if ( v == "YMD" ) globals::write_annot_date_format = YMD;
      return;
    }

  

  // read different format from EDF header
  if ( Helper::iequals( tok0 , "edf-date-format" ) )
    {
      const std::string v = Helper::toupper( tok1 );
      if ( v == "MDY" )	globals::read_edf_date_format = MDY;
      else if ( v == "DMY" ) globals::read_edf_date_format = DMY;
      else if ( v == "YMD" ) globals::read_edf_date_format = YMD;
      return;
    }

  

  // force start time/dates
  if ( Helper::iequals( tok0 , "starttime" ) )
    {
      globals::force_starttime = tok1;
      if ( globals::force_starttime.size() > 8 )
	Helper::halt( "starttime cannot be over 8 characters" );
      return;
    }
  
  
  if ( Helper::iequals( tok0 , "startdate" ) )
    {
      globals::force_startdate = tok1;
      if ( globals::force_startdate.size() > 8 )
	Helper::halt( "startdate cannot be over 8 characters" );      
      return;
    }

  
  // specify indiv-wildcard (other than ^)
  // which is needed if file ID actually contains ^
  if ( Helper::iequals( tok0 , "wildcard" ) )
    {
      globals::indiv_wildcard = tok1;
      return;
    }

  
  // always sanitize labels (channels, annots) on first reading?
  if ( Helper::iequals( tok0 , "sanitize" ) )
    {
      globals::sanitize_everything = Helper::yesno( tok1 );
      return;
    }

  
  // "auto-correct" truncated/over-long EDFs
  if ( Helper::iequals( tok0, "fix-edf" ) )
    {
      globals::autofix_edf = Helper::yesno( tok1 );
      return;
    }

  
 
  // force dig min/max (if it is invalid)
  if ( Helper::iequals( tok0 , "force-digital-minmax" ) )
    {
      globals::force_digital_minmax = Helper::yesno( tok1 );
      return;
    }

    
  // force dig min/max (if it is invalid)
  if ( Helper::iequals( tok0 , "force-digital-min" ) )
    {
      if ( ! Helper::str2int( tok1 , &globals::force_digital_min ) )
	Helper::halt( "expecting integer for force-digital-min=N" );
      return;
    }

  
    
  // force dig min/max (if it is invalid)
  if ( Helper::iequals( tok0 , "force-digital-max" ) )
    {
      if ( ! Helper::str2int( tok1 , &globals::force_digital_max ) )
	Helper::halt( "expecting integer for force-digital-max=N" );
      return;
    }
  

  
  // dp for time output
  if ( Helper::iequals( tok0, "sec-dp" ) )
    {
      if ( ! Helper::str2int( tok1 , &globals::time_format_dp ) )
	Helper::halt( "expecting integer for sec-dp=N" );
      return;
    }


  // swap spaces
  if ( Helper::iequals( tok0 , "spaces" ) )
    {
      // these are T by default, so leave as is
      // globals::replace_channel_spaces = true;
      // globals::replace_annot_spaces = true;
      if ( tok1.size() != 1 ) Helper::halt( "expecting single character after spaces" );
      globals::space_replacement = tok1[0];
      return;
    }



  // set channel names as all UPPERCASE
  if ( Helper::iequals( tok0 , "upper" ) )
    {
      globals::uppercase_channels = Helper::yesno( tok1 );
      return;
    }


    
  // if mapping to primary, retain original channel case (default=T)
  if ( Helper::iequals( tok0 , "retain-case" ) )
    {
      globals::retain_alias_case = Helper::yesno( tok1 );
      return;
    }



  // keep spaces
  if ( Helper::iequals( tok0 , "keep-spaces" ) )
    {
      globals::replace_channel_spaces = ! Helper::yesno( tok1 );
      globals::replace_annot_spaces = ! Helper::yesno( tok1 );
      return;
    }

  // keep spaces (annots only) 
  if ( Helper::iequals( tok0 , "keep-annot-spaces" ) )
    {
      globals::replace_annot_spaces = ! Helper::yesno( tok1 );
      return;
    }

      
  // keep spaces (channels only) 
  if ( Helper::iequals( tok0 , "keep-channel-spaces" ) )
    {
      globals::replace_channel_spaces = ! Helper::yesno( tok1 );
      return;
    }

  
    
  // on WRITE-ANNOTS .annot (only), set 0-duration intervals to '...' markers
  if ( Helper::iequals( tok0 , "add-ellipsis" ) ) 
    {
      globals::set_0dur_as_ellipsis =  Helper::yesno( tok1 );;
      return;
    }

  

  // treatment of gaps going from EDF+D to EDF in annots
  if ( Helper::iequals( tok0 , "annot-segment" ) )
    {
      globals::annot_disc_segment = tok1[0];
      return;
    }

      
  if ( Helper::iequals( tok0 , "annot-gap" ) )
    {
      globals::annot_disc_gap = tok1[0];
      return;
    }

  // not used??
  if ( Helper::iequals( tok0 , "annot-span-gaps" ) )
    {
      globals::annot_disc_drop_spanning = ! Helper::yesno( tok1 );
      return;
    }

  
  // split class/annot remappings (ABC/DEF|XYZ)
  if ( Helper::iequals( tok0 , "class-instance-delimiter" ) )
    {
      if ( tok1 != "" ) globals::class_inst_delimiter = tok1[0];
      return;
    }

  
    
  // combine annot class and instance IDs
  if ( Helper::iequals( tok0 , "combine-annots" ) )
    {
      globals::combine_annot_class_inst = true;
      if ( tok1 != "" ) globals::annot_class_inst_combiner = tok1[0];
      return;
    }



  // skip annots not on the whitelist (remap list)
  // also applies to EDF Annots 
  if ( Helper::iequals( tok0 , "annot-whitelist" ) )
    {
      nsrr_t::whitelist = Helper::yesno( tok1 );
      return;      
    }


  
  // skip annots on the whitelist (i.e. to get unampped remapings are a whitelist
  if ( Helper::iequals( tok0 , "annot-unmapped" ) )
    {
      nsrr_t::unmapped = Helper::yesno( tok1 );
      return;
    }
  

     
  // sleep stage prefix
  if (  Helper::iequals( tok0 , "ss-prefix" ) )
    {
      globals::sleep_stage_prefix = tok1; 
      return;
    }

  

  // POPS stages: pN1, pN2, ... 
  if ( Helper::iequals( tok0 , "ss-pops" ) )
    {
      if ( Helper::yesno( tok1 ) ) 
	globals::sleep_stage_prefix = "p";
      return;
    }

  
    
  // SOAP stages: sN1, sN2, ... 
  if ( Helper::iequals( tok0 , "ss-soap" ) )
    {
      if ( Helper::yesno( tok1 ) )
	globals::sleep_stage_prefix = "s";
      return;
    }

  

  // assume 0-dur sleep stage annots are of epoch-duration (change on loading)
  // (by default, true)
  if (  Helper::iequals( tok0 , "assume-stage-duration" ) )
    {
      globals::sleep_stage_assume_epoch_duration = Helper::yesno( tok1 );
      return;
    }

  
  
  if (  Helper::iequals( tok0 , "default-starttime" ) )
    {
      globals::use_default_starttime = true; //default anyway
      globals::default_starttime = tok1;
      return;
    }

  
  
  if (  Helper::iequals( tok0 , "no-default-starttime" ) )
    {
      globals::use_default_starttime = ! Helper::yesno( tok1 ); 
      return;
    }

  

  if (  Helper::iequals( tok0 , "default-startdate" ) )
    {
      globals::use_default_startdate = true; //default anyway
      globals::default_startdate = tok1;
      return;
    }
  
  

  if (  Helper::iequals( tok0 , "no-default-startdate" ) )
    {
      globals::use_default_startdate = ! Helper::yesno( tok1 ); 
      return;
    }

  
    
  // individual-specific variables
  if ( Helper::iequals( tok0 , "vars" ) ) 
    {
      cmd_t::attach_ivars( tok1 );
      return;
    }

  
    
  // ID re-mapper?
  if ( Helper::iequals( tok0 , "ids" ) )
    {
      cmd_t::attach_idmapper( tok1 );
      return;
    }


  
    
  // channel type labels: partial match
  if ( Helper::iequals( tok0 , "ch-match" ) )
    {
      //  type|label1|label2,type|label1|label2
      std::vector<std::string> tok2 = Helper::quoted_parse( tok1 , "," );
      for (int s=0;s<tok2.size();s++)
        {
	  std::vector<std::string> tok3 = Helper::quoted_parse( tok2[s] , "|" );
	  if ( tok3.size() < 2 ) Helper::halt( "bad format for " + tok0 + "=" + tok1 );
	  for ( int j=1;j<tok3.size();j++) globals::add_channel_map( tok3[j] , tok3[0] );
	}
      return;
    }

  
    
  // channel type labels: exact match
  if ( Helper::iequals( tok0 , "ch-exact" ) )
    {
      //  type|label1|label2,type|label1|label2
      std::vector<std::string> tok2 = Helper::quoted_parse( tok1 , "," );
      for (int s=0;s<tok2.size();s++)
        {
	  std::vector<std::string> tok3 = Helper::quoted_parse( tok2[s] , "|" );
	  if ( tok3.size() < 2 ) Helper::halt( "bad format for " + tok0 + "=" + tok1 );
	  for ( int j=1;j<tok3.size();j++) globals::add_channel_map_exact( tok3[j] , tok3[0] );
	}
      return;
    }

  

  // wipe channel type map
  if ( Helper::iequals( tok0 , "ch-clear" ) )
    {
      if ( Helper::yesno( tok1 ) ) globals::clear_channel_map();
      return;
    }

  
    
  // naughty list?
  if ( Helper::iequals( tok0 , "fail-list" ) )
    {
      globals::write_naughty_list = true;
      globals::naughty_list = tok1;
      // create an empty file (i.e. as we append things to this subsequently
      std::ofstream P( globals::naughty_list.c_str() , std::ios::out );
      P.close();      
      return;
    }

  

  // -t output compression 
  if ( Helper::iequals( tok0 , "compressed" ) )
    {
      bool yesno = Helper::yesno( tok1 );      
      globals::cmddefs().all_compressed( yesno );
      globals::cmddefs().none_compressed( !yesno );
      return;
    }


  //          default
  // stages   Y
  // others   N
  // i.e. order of annot-remap=F and nsrr-remap=T will matter

  
  
  // annot-remap: if F, then wipe all (stages + any added NSRR terms)
  if ( Helper::iequals( tok0 , "annot-remap" ) )
    {
      
      //nsrr_t::do_remap = Helper::yesno( tok1 ) ;

      // clear ALL (stages + others) pre-populated NSRR remapping      
      if ( !  Helper::yesno( tok1 ) )  
	nsrr_t::clear();
      
      return;
    }

  
  
  // nsrr-remap: if T, add in extra terms (off by default)  
  if ( Helper::iequals( tok0 , "nsrr-remap" ) )
    {
      if ( Helper::yesno( tok1 ) )
	nsrr_t::init_nsrr_mappings();
      return;
    }
  
  
  
  // generic annotation re-labelling, same format as 'alias'
  if ( Helper::iequals( tok0 , "remap" ) )
    {
      nsrr_t::annot_remapping( globals::sanitize_everything ? Helper::sanitize( tok1 ) : tok1 );
      return;
    }

  
			      
  // for EDF-annots only, set these to be a class rather than an annotation
  // (and apply any remappings)
  // if annot-whitelist=T then we *only* add EDF Annots if they are named here 
  if ( Helper::iequals( tok0 , "edf-annot-class" ) )
    {
      nsrr_t::edf_annot_class( globals::sanitize_everything ? Helper::sanitize( tok1 ) : tok1 );
      return;
    }

  
    
  if (  Helper::iequals( tok0 , "edf-annot-class-all" ) )
    {
      // equals 'edf-annot-class=*'
      if ( Helper::yesno( tok1 ) )
	nsrr_t::edf_annot_class( "*" ); // set all to be read as a class, e.g. for Moonlight
      return;
    }
  
  
  

  // fix delimiter to tab only for .annot
  // default T --> tab-only=F is option to allow spaces  
  if ( Helper::iequals( tok0 , "tab-only" ) )
    {
      globals::allow_space_delim = ! Helper::yesno( tok1 );
      return;
    }

  
  
  // if annot INST ID black, add hh:mm:ss
  if ( Helper::iequals( tok0 , "inst-hms" ) )
    {
      globals::set_annot_inst2hms = Helper::yesno( tok1 );
      return;
    }

  
    
  // set INST ID to hh:mm:ss, whether it is blank or not
  if ( Helper::iequals( tok0 , "force-inst-hms" ) )
    {
      globals::set_annot_inst2hms_force = Helper::yesno( tok1 );
      return;
    }

  
    
  // not enforce epoch check for .eannot
  // default = 5 ... (arbitrary, but allow the occassional off-by-one issue)
  if ( Helper::iequals( tok0 , "epoch-check" ) )
    {
      if ( ! Helper::str2int( tok1 , &globals::enforce_epoch_check ) )
        Helper::halt( "epoch-check requires integer value, e.g. epoch-check=10" );
      globals::enforce_epoch_check = abs( globals::enforce_epoch_check );
      return;
    }

  
    
  // set default epoch length
  if ( Helper::iequals( tok0 , "epoch-len" ) )
    {
      if ( ! Helper::str2int( tok1 , &globals::default_epoch_len ) )
	Helper::halt( "epoch-len requires integer value, e.g. epoch-len=10" );
      return;
    }

  

  // additional annot files to add from the command line
  // i.e. so we don't have to edit the sample-list
  if ( Helper::iequals( tok0 , "annot-file" ) )
    {
      globals::annot_files = Helper::parse( tok1 , "," );
      return;
    }
  
  
 
  // specified annots (only load these)
  if ( Helper::iequals( tok0 , "annot" ) ) 
    {
      param_t dummy;     
      dummy.add( "dummy" , globals::sanitize_everything ? Helper::sanitize( tok1 ) : tok1 );
      if ( tok1 == "." )
	globals::specified_annots.clear();
      else
	globals::specified_annots = Helper::combine( globals::specified_annots, dummy.strset( "dummy" , "," ) );      
      return;
    }

  
    
  // specified annots (only load these - but do not apply sanitization)
  if ( Helper::iequals( tok0 , "raw-annot" ) ) 
    {
      param_t dummy;     
      dummy.add( "dummy" , tok1 );
      if ( tok1 == "." )
	globals::specified_annots.clear();
      else
	globals::specified_annots = Helper::combine( globals::specified_annots, dummy.strset( "dummy" , "," ) );
      return;
    }
  
  // excluded annots (do not load these)
  if ( Helper::iequals( tok0 , "ignore-annot" ) ) 
    {
      param_t dummy;     
      dummy.add( "dummy" , globals::sanitize_everything ? Helper::sanitize( tok1 ) : tok1 );
      if ( tok1 == "." )
	globals::excluded_annots.clear();
      else
	globals::excluded_annots = Helper::combine( globals::excluded_annots, dummy.strset( "dummy" , "," ) );      
      return;
    }

  // specified annots (do not load these - but do not apply sanitization)
  if ( Helper::iequals( tok0 , "ignore-raw-annot" ) ) 
    {
      param_t dummy;     
      dummy.add( "dummy" , tok1 );
      if ( tok1 == "." )
	globals::excluded_annots.clear();
      else
	globals::excluded_annots = Helper::combine( globals::excluded_annots, dummy.strset( "dummy" , "," ) );
      return;
    }


    
  // delimiter char for annot key=value pairs (default '=')
  if ( Helper::iequals( tok0 , "annot-keyval" ) )
    {
      globals::annot_keyval_delim = tok1[0];
      return;
    }


  
  // alternate char for default (primary) annot meta delims (default |)
  if ( Helper::iequals( tok0 , "annot-meta-delim1" ) )
    {
      globals::annot_meta_delim = tok1[0];
      return;
    }


  // char for annot meta delims (default ; )
  if ( Helper::iequals( tok0 , "annot-meta-delim2" ) )
    {
      // i.e. setting this by default means | or this value
      globals::annot_meta_delim2 = tok1[0];
      return;
    }


  
  // default type if numeric if meta-data type of not otherwise defined
  // true by default - will give error for text, in which case need to specify  
  if ( Helper::iequals( tok0 , "annot-meta-default-num" ) )
    {
      globals::annot_default_meta_num_type = Helper::yesno( tok1 );
      return;
    }


			  
  // default type for an annot meta data (all .annots)
  if ( Helper::iequals( tok0 , "num-atype" ) )
    {
      std::vector<std::string> toka = Helper::parse( tok1 , "," );
      for (int i=0; i<toka.size(); i++)
	globals::atypes[ toka[i] ] = globals::A_DBL_T;
      return;
    }


    
  if ( Helper::iequals( tok0 , "txt-atype" ) || Helper::iequals( tok0 , "str-atype" ) )
    {
      std::vector<std::string> toka = Helper::parse( tok1 , "," );
      for (int i=0; i<toka.size(); i++)
	globals::atypes[ toka[i] ] = globals::A_TXT_T;
      return;
    }
  


  if ( Helper::iequals( tok0 , "int-atype" ) )
    {
      std::vector<std::string> toka = Helper::parse( tok1 , "," );
      for (int i=0; i<toka.size(); i++)
	globals::atypes[ toka[i] ] = globals::A_INT_T;
      return;
    }


  
  if ( Helper::iequals( tok0 , "bool-atype" ) )
    {
      std::vector<std::string> toka = Helper::parse( tok1 , "," );
      for (int i=0; i<toka.size(); i++)
	globals::atypes[ toka[i] ] = globals::A_BOOL_T;
      return;
    }


  
  // annotation alignment
  if ( Helper::iequals( tok0 , "align-annots" ) )
    {
      globals::annot_alignment = Helper::vec2set( Helper::parse( tok1 , "," ) ) ;
      return;
    }
  

    
  // signal alias?
  if ( Helper::iequals( tok0 , "alias" ) )
    {

      const std::string str = globals::sanitize_everything ?
	(  globals::replace_channel_spaces
	   ? Helper::trim( Helper::sanitize( tok1 ) , '_' ) 
	   : Helper::trim( Helper::sanitize( tok1 , ' ' ) , '_' ) 
	   )
	: tok1 ; 
      cmd_t::signal_alias( str );
      return;
    }


  
  // behaviour when a problem encountered
  if ( Helper::iequals( tok0 , "bail-on-fail" ) )
    {
      globals::bail_on_fail = Helper::yesno( tok1 );
      return;
    }



  // force reading EDF+ as EDF
  // (skips EDF Annotations track, and ignores any time-track information)
  if ( Helper::iequals( tok0 , "force-edf" ) )
    {
      globals::force_edf = Helper::yesno( tok1 );      
      return;
    }

  

  
  // pre-read whole EDF(+D) rather than random access reads
  //  (to avoid slow fseek() implementations) 
  if ( Helper::iequals( tok0, "preload" ) ) 
    {
      globals::edf_stream_read = Helper::yesno( tok1 );
      return;
    }
  

  
  // skip anyt EDF Annotations from EDF+
  if ( Helper::iequals( tok0 , "skip-edf-annots" ) )
    {
      globals::skip_edf_annots = Helper::yesno( tok1 );
      return;
    }



  // do not load sample-list annotations
  if ( Helper::iequals( tok0 , "skip-sl-annots" ) )
    {
      globals::skip_sl_annots = Helper::yesno( tok1 );
      return;
    }


    
  // do not read ANY annotations
  if ( Helper::iequals( tok0 , "skip-annots" ) ||
       Helper::iequals( tok0 , "skip-all-annots" ) )
    {
      // nb - internally skip_edf_annots and skip_sl_annots are redundant I believe, but can keep as is 
      globals::skip_edf_annots = globals::skip_sl_annots = globals::skip_nonedf_annots = Helper::yesno( tok1 );
      return;
    }


    
  // project path
  if ( Helper::iequals( tok0 , "path" ) )
    {
      globals::param.add( "path" , tok1 );
      return;
    }



  // prepend/append for text-table output files
  if ( Helper::iequals( tok0 , "tt-prepend" ) ||  Helper::iequals( tok0 , "tt-prefix" ) )
    {
      globals::txt_table_prepend = tok1;
      return;
    }



  if ( Helper::iequals( tok0 , "tt-append" ) ||  Helper::iequals( tok0 , "tt-suffix" ) )
    {
      globals::txt_table_append = tok1;
      return;
    }


  
  // shift times +12 hours if between 04:00 and 12:00  (--> 16:00 to 00:00 ) 
  if ( Helper::iequals( tok0 , "assume-pm-start" ) )
    {

      if ( tok1 == "0"
	   || Helper::iequals( tok1 , "n" )
	   || Helper::iequals( tok1 , "no" ) )
	globals::assume_pm_starttime = false; 
      else
	{
	  globals::assume_pm_starttime = true;
	  
	  if ( ! Helper::str2int( tok1 , &globals::assume_pm_starttime_hour ) ) 
	    Helper::halt( "expecting integer between 0 and 12" );

	}
      
      return;
    }

  // power band defintions
  if ( Helper::iequals( tok0 , "slow" ) 
       || Helper::iequals( tok0 , "delta" ) 
       || Helper::iequals( tok0 , "theta" ) 
       || Helper::iequals( tok0 , "alpha" ) 
       || Helper::iequals( tok0 , "sigma" ) 
       || Helper::iequals( tok0 , "beta" ) 
       || Helper::iequals( tok0 , "gamma" ) 
       || Helper::iequals( tok0 , "total" ) ) 
    {
      std::vector<std::string> f = Helper::parse( tok1 , ",-" );
      if ( f.size() != 2 ) Helper::halt( "expecting band=lower,upper" );
      double f0, f1;
      if ( ! Helper::str2dbl( f[0] , &f0 ) ) Helper::halt( "expecting numeric for power range" );
      if ( ! Helper::str2dbl( f[1] , &f1 ) ) Helper::halt( "expecting numeric for power range" );
      if ( f0 >= f1 ) Helper::halt( "expecting band=lower,upper" );
      if ( f0 < 0 || f1 < 0 ) Helper::halt( "negative frequencies specified" );
      
      if      ( Helper::iequals( tok0 , "slow" ) )  globals::freq_band[ SLOW ] =  freq_range_t( f0 , f1 ) ;
      else if ( Helper::iequals( tok0 , "delta" ) ) globals::freq_band[ DELTA ] = freq_range_t( f0 , f1 ) ;
      else if ( Helper::iequals( tok0 , "theta" ) ) globals::freq_band[ THETA ] = freq_range_t( f0 , f1 ) ;
      else if ( Helper::iequals( tok0 , "alpha" ) ) globals::freq_band[ ALPHA ] = freq_range_t( f0 , f1 ) ;
      else if ( Helper::iequals( tok0 , "sigma" ) ) globals::freq_band[ SIGMA ] = freq_range_t( f0 , f1 ) ;
      else if ( Helper::iequals( tok0 , "beta" ) )  globals::freq_band[ BETA ] =  freq_range_t( f0 , f1 ) ;
      else if ( Helper::iequals( tok0 , "gamma" ) ) globals::freq_band[ GAMMA ] = freq_range_t( f0 , f1 ) ;
      else if ( Helper::iequals( tok0 , "total" ) ) globals::freq_band[ TOTAL ] = freq_range_t( f0 , f1 ) ;
      
    }


  // midflight ... allow 10 user-defined bands, but these are not
  // yet plugged into bandaid_t... TODO
  if ( Helper::iequals( tok0 , "band" ) )
    {
      // expect band=BAND1,lwr,upr
      std::vector<std::string> tok = Helper::parse( tok1 , ",-" );
      if ( tok.size() != 3 ) Helper::halt( "expecting band=label,lwr,upr" );
      std::string bstr = Helper::toupper( tok[0] );
      double lwr, upr;
      if ( ! Helper::str2dbl( tok[1] , &lwr ) ) Helper::halt( "lwr is not numeric" );
      if ( ! Helper::str2dbl( tok[2] , &upr ) ) Helper::halt( "upr is not numeric" );
      if ( lwr >= upr ) Helper::halt( "lwr >= upr" );
      if ( lwr < 0 ) Helper::halt( "frequencies must be positive" );
      frequency_band_t b = globals::band( bstr );
      if ( b == UNKNOWN_BAND ) Helper::halt( "band label not recognized" );
      globals::set_band( b , lwr , upr );
      return;
    }


  
  // exclude individuals?
  if ( Helper::iequals( tok0 , "exclude" ) )
    {
      
      if ( globals::id_includes.size() > 0 ) 
	Helper::halt( "cannot specify both include= and exclude= lists" );

      std::string xfile = Helper::expand( tok1 );
      
      if ( Helper::fileExists( xfile ) ) 
	{
	  std::ifstream XIN( xfile.c_str() , std::ios::in );
	  while ( !XIN.eof() ) 
	    {
	      // format: ID {white-space} any notes (ignored)
	      std::string line2;
	      Helper::safe_getline( XIN , line2);
	      if ( XIN.eof() || line2 == "" ) continue;
	      std::vector<std::string> tok2 = Helper::parse( line2 , "\t " );
	      if ( tok2.size() == 0 ) continue;			      
	      std::string xid = tok2[0];
	      globals::id_excludes.insert( xid );
	    }

	  logger << "  excluding " << globals::id_excludes.size() 
		 << " individuals from " << xfile << "\n";
	  XIN.close();
	}
      else 
	Helper::halt( "exclude file " + xfile + " does not exist" );

      return;
    }


  
  // exclude individuals?
  if ( Helper::iequals( tok0 , "include" ) )
    {
      
      if ( globals::id_excludes.size() > 0 ) 
	Helper::halt( "cannot specify both include= and exclude= lists" );

      std::string xfile = Helper::expand( tok1 );
      
      if ( Helper::fileExists( xfile ) ) 
	{
	  std::ifstream XIN( xfile.c_str() , std::ios::in );
	  while ( !XIN.eof() ) 
	    {
	      // format: ID {white-space} any notes (ignored)
	      std::string line2;
	      Helper::safe_getline( XIN , line2);
	      if ( XIN.eof() || line2 == "" ) continue;
	      std::vector<std::string> tok2 = Helper::parse( line2 , "\t " );
	      if ( tok2.size() == 0 ) continue;			      
	      std::string xid = tok2[0];
	      globals::id_includes.insert( xid );
	    }

	  logger << "  only including " << globals::id_includes.size() 
		 << " individuals from " << xfile << "\n";
	  XIN.close();
	}
      else 
      	Helper::halt( "include file " + xfile + " does not exist" );

      return;
    }


  // not sure if/where this is used now
  
  if ( tok0[0] == '-' )
    {
      globals::param.add( tok0.substr(1) , tok1 );
      return;
    }

  

  
}


//
// Channel type information
//

void cmd_t::define_channel_type_variables( edf_t & edf )
{
  // add these even if blank, so that scripts can always use
  
  std::string eeg = globals::list_channels( EEG , edf.header.label );
  cmd_t::ivars[ edf.id ][ "eeg" ] = eeg;
    
  std::string ref = globals::list_channels( REF , edf.header.label );
  cmd_t::ivars[ edf.id ][ "ref" ] = ref;

  std::string ic = globals::list_channels( IC , edf.header.label );
  cmd_t::ivars[ edf.id ][ "ic" ] = ic;

  std::string imf = globals::list_channels( IMF , edf.header.label );
  cmd_t::ivars[ edf.id ][ "imf" ] = imf;

  std::string eog = globals::list_channels( EOG , edf.header.label );
  cmd_t::ivars[ edf.id ][ "eog" ] = eog;

  std::string ecg = globals::list_channels( ECG , edf.header.label );
  cmd_t::ivars[ edf.id ][ "ecg" ] = ecg;

  std::string emg = globals::list_channels( EMG , edf.header.label );
  cmd_t::ivars[ edf.id ][ "emg" ] = emg;

  std::string leg = globals::list_channels( LEG , edf.header.label );
  cmd_t::ivars[ edf.id ][ "leg" ] = leg;

  std::string generic = globals::list_channels( GENERIC , edf.header.label );
  cmd_t::ivars[ edf.id ][ "generic" ] = generic;
  
  std::string airflow = globals::list_channels( AIRFLOW , edf.header.label );
  cmd_t::ivars[ edf.id ][ "airflow" ] = airflow;
  
  std::string effort = globals::list_channels( EFFORT , edf.header.label );
  cmd_t::ivars[ edf.id ][ "effort" ] = effort;

  std::string oxygen = globals::list_channels( OXYGEN , edf.header.label );
  cmd_t::ivars[ edf.id ][ "oxygen" ] = oxygen;

  std::string position = globals::list_channels( POSITION , edf.header.label );
  cmd_t::ivars[ edf.id ][ "position" ] = position;
  
  std::string light = globals::list_channels( LIGHT , edf.header.label );
  cmd_t::ivars[ edf.id ][ "light" ] = light;
  
  std::string snore = globals::list_channels( SNORE , edf.header.label );
  cmd_t::ivars[ edf.id ][ "snore" ] = snore;
  
  std::string hr = globals::list_channels( HR , edf.header.label );
  cmd_t::ivars[ edf.id ][ "hr" ] = hr;
  
  std::string ignore = globals::list_channels( IGNORE_SIGNAL , edf.header.label );
  cmd_t::ivars[ edf.id ][ "ignore" ] = ignore;

  
}


//
// Kludge... so that scope compiles on Mac... need to ensure the reduce_t is 
// called within the primary libluna.dylib, or else this function is not accessible
// when linking against libluna.dylib from scope...   TODO... remind myself 
// how to control which functions are exported, etc, in shared libraries
//

void tmp_includes()
{

  // dummy function to make sure reduce_t() [ which otherwise isn't used by luna-base ] is 
  // accessible from scope

  std::vector<double> d;
  std::vector<uint64_t> tp;
  uint64_t s1,s2;
  reduce_t r( &d,&tp,s1,s2,1);
  return;
}



//
// Make an EDF continuous
// 

void proc_continuous( edf_t & edf , param_t & param )
{
  logger << " forcing EDF to be continuous\n";
  edf.set_edf();
}


//
// Output DB 
//

std::string cmd_t::resolved_outdb( const std::string & id , const std::string & str )
{ 
  return Helper::insert_indiv_id( Helper::sanitize( id ) , str ); 
} 



//
// Attach i-vars from one or more files
//

void cmd_t::attach_ivars( const std::string & file )
{
  
  std::vector<std::string> files = Helper::parse( file , "," );
  
  for ( int f = 0; f < files.size() ; f++ )
    {
      std::string filename = Helper::expand( files[f] );

      if ( ! Helper::fileExists( filename ) )
	Helper::halt( "could not find " + filename );
      
      //      logger << "reading from " << filename << "\n";
      
      std::ifstream IN1( filename.c_str() , std::ios::in );
      bool header = true ; 
      int idcol = -1;
      int ncols = 0;

      std::vector<std::string> head;

      while ( ! IN1.eof() ) 
	{
          std::string s;
          Helper::safe_getline( IN1 , s );
          if ( IN1.eof() ) break;
          if ( s == "" ) continue;

	  std::vector<std::string> tok = Helper::parse( s , "\t" );
	  
	  if ( header )
	    {
	      for (int i=0;i<tok.size();i++)
		{
		  if ( tok[i] == "ID" )
		    {
		      if ( idcol != -1 )
			Helper::halt( "cannot have multiple ID columns in " + filename );
		      idcol = i;
		    }
		}
	      // store header
	      head = tok;
	      ncols = head.size();	      
	    }
	  
	  if ( ! header )
	    {
	      if ( ncols != tok.size() )
		{
		  std::cerr << " expecting " << ncols << " columns from header\n";
		  std::cerr << " observed " << tok.size() << "\n"
			    << s << "\n";
		  Helper::halt( "inconsistent number of columns in " + filename );
		}
	      for (int c=0;c<ncols;c++)
		{
		  if ( c == idcol ) continue;		  
		  cmd_t::ivars[ tok[ idcol ] ][ head[ c ] ] = tok[c] ;		  
		  //logger << "setting: " << tok[ idcol ] << " " << head[ c ]  << "  " << tok[c]  << "\n";
		  
		}
	      
	    }

	  // read header now
	  header = false; 
	}
      IN1.close();
      
      //logger << "  attached " << ncols - 1 << " from " << filename << "\n";
    }
}



//
// Attach idmapper from a file
//

void cmd_t::attach_idmapper( const std::string & file )
{

  std::string filename = Helper::expand( file );

  if ( ! Helper::fileExists( filename ) )
    Helper::halt( "could not find " + filename );
      
  std::ifstream IN1( filename.c_str() , std::ios::in );
  while ( ! IN1.eof() ) 
    {
      std::string s;
      Helper::safe_getline( IN1 , s );
      if ( IN1.eof() ) break;
      if ( s == "" ) continue;
      std::vector<std::string> tok = Helper::parse( s , "\t" );
      if ( tok.size() != 2 ) Helper::halt( "bad format in " + filename );
      cmd_t::idmapper[ tok[0] ] = tok[1];
    }

  IN1.close();

  logger << "  read " << cmd_t::idmapper.size() << " IDs to remap\n";

}


void cmd_t::register_specials()
{

  specials.insert( "ch-match" );
  specials.insert( "ch-exact" ) ;
  specials.insert( "ch-clear" ) ;

  // 14 channel types:: automatic, but okay to overwrite...  
  // e.g. $eeg $emg $leg etc 

  specials.insert( "spaces" );
  specials.insert( "keep-spaces" );
  specials.insert( "keep-annot-spaces" );
  specials.insert( "keep-channel-spaces" );  
  specials.insert( "silent" ) ;
  specials.insert( "id" );
  specials.insert( "verbose" ) ;
  specials.insert( "devel" );
  specials.insert( "sec-dp" );
  specials.insert( "sig" ) ;
  specials.insert( "vars" );
  specials.insert( "ids" );    
  specials.insert( "add" ) ;
  specials.insert( "ss-prefix" );

  specials.insert( "fail-list" ) ;
  specials.insert( "compressed" ) ;
  specials.insert( "nsrr-remap" ) ;
  specials.insert( "remap" ) ;
  specials.insert( "combine-annots");
  specials.insert( "class-instance-delimiter");
  specials.insert( "tab-only" );
  specials.insert( "annot-folder" ) ;
  specials.insert( "annots-folder" ) ; 
  specials.insert( "inst-hms" ) ;
  specials.insert( "force-inst-hms" ) ;
  specials.insert( "no-epoch-check" ) ;
  specials.insert( "epoch-len" ) ;
  specials.insert( "annots-file" ) ;
  specials.insert( "annots-files" ) ;
  specials.insert( "annot-file" ) ;
  specials.insert( "annot-files" ) ;
  specials.insert( "annots" );
  specials.insert( "annot"  );
  specials.insert( "alias" ) ;
  specials.insert( "bail-on-fail" ) ;
  specials.insert( "force-edf" ) ;
  specials.insert( "skip-edf-annots" ) ;
  specials.insert( "skip-annots" ) ;
  specials.insert( "skip-all-annots" ) ;
  specials.insert( "path" ) ;
  specials.insert( "tt-prepend" );
  specials.insert( "tt-prefix" ) ;
  specials.insert( "tt-append" ); 
  specials.insert( "tt-suffix" );
  specials.insert( "assume-pm-start" );
  specials.insert( "slow" ) ;
  specials.insert( "delta" ) ;
  specials.insert( "theta" ) ;
  specials.insert( "alpha" ) ;
  specials.insert( "sigma" ) ;
  specials.insert( "beta" ) ;
  specials.insert( "gamma" ) ;
  specials.insert( "total" ) ;
  specials.insert( "exclude" ) ;
  specials.insert( "include" ) ;


  // register and define these topographical special variables:
  // specials.insert( "left" );
  // specials.insert( "midline" );
  // specials.insert( "right" );
  // specials.insert( "anterior" );
  // specials.insert( "central" );
  // specials.insert( "posterior" );
  // specials.insert( "anterio-frontal" );  
  // specials.insert( "mid-central" );  
  // specials.insert( "centro-parietal" );  
  // specials.insert( "frontal" );  
  // specials.insert( "fronto-central" ); 
  // specials.insert( "occiptital" );  
  // specials.insert( "parietal" );  
  // specials.insert( "parieto-occipital" );  
  // specials.insert( "pre-frontal" );  
  // specials.insert( "temporal" );

  //
  // EEG topographical groupings
  //

  //
  // Left/right
  //

  vars[ "left" ] 
    = "FP1,AF7,AF3,F1,F3,F5,F7,FT7,FC5,FC3,FC1,C1,C3,C5,T7,TP7,CP5,CP3,CP1,P1,P3,P5,P7,P9,PO7,PO3,O1";
  
  vars[ "midline" ] 
    = "IZ,OZ,POZ,PZ,CPZ,FPZ,AFZ,FZ,FCZ,CZ";
  
  vars[ "right" ] 
    = "FP2,AF8,AF4,F2,F4,F6,F8,FT8,FC6,FC4,FC2,C2,C4,C6,T8,TP8,CP6,CP4,CP2,P2,P4,P6,P8,P10,PO8,PO4,O2";

  // 
  // Anterior/posterior
  //

  vars[ "anterior" ] 
    = "FP1,AF7,AF3,F1,F3,F5,F7,FPZ,AFZ,FZ,FP2,AF8,AF4,F2,F4,F6,F8";

  vars[ "central" ] 
    = "FT7,FC5,FC3,FC1,C1,C3,C5,T7,TP7,CP5,CP3,CP1,CPZ,FCZ,CZ,FT8,FC6,FC4,FC2,C2,C4,C6,T8,TP8,CP6,CP4,CP2";
  
  vars[ "posterior" ] 
    = "P1,P3,P5,P7,P9,PO7,PO3,O1,IZ,OZ,POZ,PZ,P2,P4,P6,P8,P10,PO8,PO4,O2";
  

  // 
  // Regions
  //

  vars[ "pre-frontal" ]  
    = "FP1,FPZ,FP2";

  vars[ "anterio-frontal" ]
    = "AF7,AF3,AFZ,AF8,AF4";

  vars[ "mid-central" ]  
    = "C1,C3,C5,CZ,C2,C4,C6";
    
  vars[ "centro-parietal" ]  
    = "CP5,CP3,CP1,CPZ,CP6,CP4,CP2";

  vars[ "frontal" ]  
    = "F1,F3,F5,F7,FZ,F2,F4,F6,F8";
  
  vars[ "fronto-central" ]  
    = "FC5,FC3,FC1,FCZ,FC6,FC4,FC2";

  vars[ "occiptital" ]  
    = "O1,IZ,OZ,O2";

  vars[ "parietal" ]  
    = "P1,P3,P5,P7,P9,PZ,P2,P4,P6,P8,P10";
  
  vars[ "parieto-occipital" ]  
    = "PO7,PO3,POZ,PO8,PO4";
    
  vars[ "temporal" ]  
    = "FT7,T7,TP7,FT8,T8,TP8";


  // specials.insert( "eeg" );
  // specials.insert( "ref" );
  // specials.insert( "ic" );
  // specials.insert( "eog" );
  // specials.insert( "ecg" );
  // specials.insert( "emg" );
  // specials.insert( "leg" );  
  // specials.insert( "generic" );
  // specials.insert( "airflow" );
  // specials.insert( "effort" );  
  // specials.insert( "oxygen" );
  // specials.insert( "position" );
  // specials.insert( "light" );
  // specials.insert( "snore" );
  // specials.insert( "hr" );
  // specials.insert( "ignore" );

}


std::map<std::string,int> cmd_t::pull_ivar( const std::vector<std::string> & ids , const std::string & phe )
{
  
  std::map<std::string,int> retval;
  // ivars:  ID -> PHE -> VAL

  for (int i=0;i<ids.size();i++)
    {
      if ( ivars.find( ids[i] ) == ivars.end() ) continue;

      const std::map<std::string,std::string> & data = ivars.find( ids[i] )->second;

      if ( data.find( phe ) == data.end() ) continue;

      // attach only valid integer values here

      int x;
      if ( Helper::str2int( data.find( phe )->second , &x ) )
	retval[ ids[i] ] = x;
    }

  return retval;
}


bool cmd_t::pull_ivar( const std::string & id , const std::string & phe , double * x )
{
  
  // ivars:  ID -> PHE -> VAL

  // ID
  if ( ivars.find( id ) == ivars.end() ) return false;

  // var
  const std::map<std::string,std::string> & data = ivars.find( id )->second;
  if ( data.find( phe ) == data.end() ) return false;

  // attach only valid doubles
  return Helper::str2dbl( data.find( phe )->second , x ); 

}

bool cmd_t::pull_ivar_bool( const std::string & id , const std::string & phe )
{
  
  // ivars:  ID -> PHE -> VAL

  // ID
  if ( ivars.find( id ) == ivars.end() ) return false;
  
  // var
  const std::map<std::string,std::string> & data = ivars.find( id )->second;
  if ( data.find( phe ) == data.end() ) return false;
  
  // return value
  return Helper::yesno( data.find( phe )->second );

}

// REQUIRES

void proc_requires( edf_t & edf , param_t & param )
{
  // luna version?
  if ( param.has( "version" ) )
    {
      //  version = "v1.2.2 ";

      std::string v = param.value( "version" );
      
      if ( param.empty( "version" ) || v.size() == 0 ) Helper::halt( "no valid version label specified" );

      if ( v[0] == 'v' || v[0] == 'V' ) v = v.substr(1);

      std::vector<std::string> tok = Helper::parse( v, "." );
      if ( tok.size() < 1 || tok.size() > 3 ) Helper::halt( "invalid version code" );

      const bool has2 = tok.size() >= 2;
      const bool has3 = tok.size() >= 3;
      
      int n1 = -1 , n2 = -1 , n3 = -1;
      
      if ( ! Helper::str2int( tok[0] , &n1 ) )
	Helper::halt( "invalid version code" );

      if ( has2 && ! Helper::str2int( tok[1] , &n2 ) )
	Helper::halt( "invalid version code" );
      
      if ( has3 && ! Helper::str2int( tok[2] , &n3 ) )
	Helper::halt( "invalid version code" );

      if ( n1 > globals::major_version_number ||
	   ( has2 && n2 > globals::minor_version_number ) ||
	   ( has3 && n3 > globals::patch_version_number ) )
	Helper::halt( "required version " + v + " but current Luna is " + globals::version );

      logger << "  current Luna version " << globals::version << " satisfies required version " << v << "\n";
      
    }
    
}

// CONTAINS

void proc_has_signals( edf_t & edf , param_t & param )
{

  // operates in different modes:  
  //    - alters return code (default)
  //    - skips to next EDF (if skip option is given) 
  //    - or set a variable

  const bool skip = param.has( "skip" ) || param.has( "skip-if-none" ) ;
  
  const bool skip_if_none = param.has( "skip-if-none" );
   
  // check this EDF has the signals OR annots OR stage information
  
  const bool check_stages = param.has( "stages" );

  const bool check_annots = param.has( "annots" ) || param.has( "annot" );

  // only signals has skip vs skip-if-none distinction
  if ( ( check_stages || check_annots ) && skip_if_none )
    Helper::halt( "cannot specify stages/annots and skip-if-none - use 'skip' instead" );
    
  
  // will always have a default signal value: if not specified, means everything
  const bool check_signals = param.value( "sig" ) != "*" ;

  if ( ( check_stages && check_annots ) ||
       ( check_stages && check_signals ) ||
       ( check_annots && check_signals ) ) 
    Helper::halt( "can only only specify stages OR annots OR sig for CONTAINS" );

  // return codes: 0   all EDFs have all signals/annots
  //               1   all EDFs have at least one of these signals/annots
  //               2   at least some EDFs do not have any of these signals/annots
  
  // 0  ALL channels in all EDFs
  // 1  only SOME channels in one or more EDFs
  // 2  NO channels in one or more EDFs
  
  // retcode starts at 0
  // if EDF has all  : do not change retcode == 0 
  // if EDF has some : retcode == 1
  // if EDF has none : retcode == 2
  // & retcode never gets smaller
  

  if ( check_stages )
    {
      
      //
      // try to make and extract stages
      // by default, this sets sslabel to SleepStage
      //
      
      edf.annotations->make_sleep_stage( edf.timeline );
      
      annot_t * annot = (*edf.annotations)( "SleepStage" );
      
      bool present = true;
      
      if ( annot == NULL )
	{

	  present = false;

	  if ( skip ) 
	    {
	      globals::problem = true;
	      return;
	    }
	  
	  // otherwise, flag a problem 
	  globals::retcode = 2;
	 	 
	}


      // Stages present, but not the correct amount?

      if ( annot != NULL )
	{
	  
	  bool has_stages = edf.timeline.hypnogram.construct( &edf.timeline , param , false ) ;

	  if ( ! has_stages ) present = false;
	  
	  if ( has_stages ) 
	    {
	      int ne = edf.timeline.num_epochs();

	      if ( ne != edf.timeline.hypnogram.stages.size() )
		{

		  present = false; 
		  
		  if ( skip ) 
		    {
		      globals::problem = true;
		      return;
		    }		  
		  
		  // otherwise, flag a problem: partial staging means 
		  if ( globals::retcode == 0 ) 
		    globals::retcode = 1;
		  
		}

	      int s_n1 = 0 , s_n2 = 0 , s_n3 = 0 , s_rem = 0 , s_wake = 0 , s_other = 0;
	      
	      const int ss = edf.timeline.hypnogram.stages.size();
	      
	      for (int s=0; s<ss; s++)
		{	  
		  
		  if      ( edf.timeline.hypnogram.stages[ s ] == WAKE ) ++s_wake;
		  else if ( edf.timeline.hypnogram.stages[ s ] == NREM1 ) ++s_n1;
		  else if ( edf.timeline.hypnogram.stages[ s ] == NREM2 ) ++s_n2;
		  else if ( edf.timeline.hypnogram.stages[ s ] == NREM3 ) ++s_n3;
		  else if ( edf.timeline.hypnogram.stages[ s ] == NREM4 ) ++s_n3;
		  else if ( edf.timeline.hypnogram.stages[ s ] == REM ) ++s_rem;
		  else ++s_other;
		}
	      
	      std::stringstream sss ;
	      sss << "N1:" << s_n1 << ","
		  << "N2:" << s_n2 << ","
		  << "N3:" << s_n3 << ","
		  << "R:" << s_rem << ","
		  << "W:" << s_wake << ","
		  << "?:" << s_other ;
	      
	      writer.value( "STAGE_COUNTS" , sss.str() );
	      
	      bool has_nrem =  ( s_n1 + s_n2 + s_n3 ) > 0 ;
	      bool has_rem  =  s_rem > 0 ;
	      bool has_wake = s_wake > 0 ;
	      
	      const int n_stages = has_rem + has_nrem + has_wake;
	      writer.value( "UNIQ_STAGES" , n_stages );
	      
	    }
	}
      
      writer.value( "STAGES" , present );
      
      return;
    }
  
  
  //
  // annotations
  //

  if ( check_annots ) 
    {

      std::vector<std::string> annots = param.has( "annot" ) ? 
	param.strvector_xsigs( "annot" ) : param.strvector_xsigs( "annots" );

      const int na = annots.size();

      int count = 0 ;

      for (int a=0; a<na; a++)
	{
	  
	  annot_t * annot = (*edf.annotations)( annots[a] ) ;
	  bool found = annot != NULL ; 	  
	  writer.level( annots[a] , globals::annot_strat );	  
	  writer.value( "PRESENT" , found );
	  if ( found ) ++count;
	}
      writer.unlevel( globals::annot_strat );

      writer.value( "NA_REQ" , na );
      writer.value( "NA_OBS" , count );
      
      // behavior based on absence?


      //
      // setting a variable?, or adjust return code
      //
      
      if ( param.has( "var" ) )
	{
	  // T = all
	  // F = ! all

	  std::string var = param.value( "var" );
	  cmd_t::ivars[ edf.id ][ var ] = count == na ? "T" : "F" ;
	  logger << "  setting ${" << var << "} <- " << ( count == na ? "T" : "F" ) << "\n";
	}

      
      //
      // in skip mode, bail if none present
      //

      if ( skip && count == 0 )
	{
	  globals::problem = true;
	  return;
	}
      
      //
      // otherwise adjust return code
      //
      
      if ( count == 0 ) globals::retcode = 2;
      else if ( count < na && globals::retcode == 0 ) globals::retcode = 1;
      
      return;
    }
  

  //
  // otherwise, check signals
  //
  

  int count = 0;

  std::vector<std::string> signals = param.strvector( "sig" );

  const int ns = signals.size();
  
  for (int s=0; s<ns; s++)
    {
      writer.level( signals[s] , globals::signal_strat );
      bool found = edf.header.has_signal( signals[s] );
      writer.value( "PRESENT" , found ); 
      if ( found ) ++count;
    }
  writer.unlevel( globals::signal_strat );

  writer.value( "NS_REQ" , ns );
  writer.value( "NS_OBS" , count );
  writer.value( "NS_TOT" , edf.header.ns );


  //
  // in skip mode
  //
  
  if ( skip )
    {
      
      int code = 0;                    // all
      if      ( count == 0 ) code = 2; // none
      else if ( count < ns ) code = 1; // some, but not all
      
      // under skip-if-none, do we have /none/?
      if ( skip_if_none && code == 2 )
	{
	  globals::problem = true;
	  return;
	}

      // else, more stringent: skip unless we have all
      if ( ( ! skip_if_none ) && code != 0 )
	{
	  globals::problem = true;
          return;
	}
      
    }


  //
  // setting a variable?, or adjust return code
  //

  if ( param.has( "var" ) )
    {
      // T = all
      // F = ! all
      std::string var = param.value( "var" );      
      cmd_t::ivars[ edf.id ][ var ] = count == ns ? "T" : "F" ;
      logger << "  setting ${" << var << "} <- " << ( count == ns ? "T" : "F" ) << "\n";
    }
  else
    {
      
      //
      // otherwise adjust return code
      //
      
      if ( count == 0 ) globals::retcode = 2;
      else if ( count < ns && globals::retcode == 0 ) globals::retcode = 1;
    } 
}

std::map<std::string,std::string>  cmd_t::indiv_var_map( const std::string & id )
{

  std::map<std::string,std::string> allvars = vars;

  if ( ivars.find( id ) != ivars.end() )
    {
      const std::map<std::string,std::string> & newvars = ivars.find( id )->second;
      std::map<std::string,std::string>::const_iterator vv = newvars.begin();
      while ( vv != newvars.end() )
	{
	  allvars[ vv->first ] = vv->second;
	  ++vv;
	}
    }
  return allvars;
}



std::vector<std::string> cmd_t::proc_all( const std::vector<std::string> & lines ,
					  std::map<std::string,std::string> * allvars )
{

  // per line:

  // swap in new vars (which may be defined on-the-fly) 
  //       Helper::swap_in_variables( &currline , &allvars );

  // expand numerics [x][y]
  //       Helper::expand_numerics( &currline );
  
    
  std::vector<std::string> r;
  
  //
  // initial checks for balanced LOOP/END-LOOP pairs in the unexpanded case
  //  -- i.e. this should always be balanced;
  //  -- we may need to check on-the-fly too, e.g. as a conditional may create issues
  //  --
  //    FOR index=a vals=1,2,3
  //
  //    IF x
  //       ...
  //       END-LOOP
  //    FI
  //
  //    IFNOT x
  //      ...
  //      END-LOOP
  //    FI
  //
  //    So best to check on-the-fly only

  std::vector<int> fdep, idep;
  
  for (int i=0; i<lines.size(); i++)
    {
      std::vector<std::string> ctok = Helper::quoted_parse( lines[i] , "\t " );
      if ( ctok.size() == 0 ) continue;
      const std::string str = Helper::toupper( ctok[0] );      

      const bool is_for = str == "LOOP";
      const bool is_until = str == "END-LOOP";
      const bool is_if = str == "IF" || str == "IFNOT";
      const bool is_fi =  str == "FI" || str == "ENDIF";

      // for a new FOR, track the depth of IF      
      if ( is_for ) fdep.push_back( idep.size() );

      // and vice versa
      if ( is_if ) idep.push_back( fdep.size() );

      // loop ending out, ensure the current opposite depth == the opening one
      if ( is_until )
	{
	  if ( fdep.size() == 0 ) Helper::halt( "unbalanced LOOP/END-LOOP pairs" );
	  if ( fdep[ fdep.size() - 1 ] != idep.size() )
	    Helper::halt( "non-nested IF/LOOP pairings" );	  
	  fdep.pop_back();
	}

      // and vice versa
      if ( is_fi )
	{
	  if ( idep.size() == 0 ) Helper::halt( "unbalanced IF/FI pairs" );
	  if ( idep[ idep.size() - 1 ] != fdep.size() )
	    Helper::halt( "non-nested IF/LOOP pairings" );	  
	  idep.pop_back();
	}      	      
    }


  //
  // count loops and ifs
  //
  
  int loop_count = 0;
  
  int if_count = 0;
  bool excise = false;
  std::string if_condition = "";
  
    
  // make expansions
  
  std::vector<std::string> keys;    // for nested loops, last is the inner
  std::map<std::string,int> slots;  // for a given loop, which slot we at?
  std::map<std::string,std::vector<std::string> > seqs; // for a given loop, what are the slots?
  std::map<std::string,std::string> vals; // swap in this value
  std::map<std::string,int> repeat_pointer; // line to loop back to 

  std::string inner = ""; // inner loop


  if ( globals::mirror )
    logger << "\n ------------------------------------------------------------\n"
	   << " mirrored command inputs:\n\n";
  
  
  //
  // Primary iteration through 
  //

  const int nl = lines.size();
  
  for (int i=0; i<nl; i++)
    {
      
      
      // line to process
      std::string l1 = lines[i];

      //
      // only do these substitutions/expansions, etc, if not
      // excised (to avoid ${var} assignments)
      //
      
      if ( ! excise ) 
	{
	  
	  //
	  // Also, swap any loop indexes (e.g. as they may feature in a LOOP vals=arg)
	  //
	  
	  std::map<std::string,std::string>::const_iterator vv = vals.begin();
	  while ( vv != vals.end() )
	    {
	      l1 = Helper::search_replace( l1 , vv->first , vv->second ) ;
	      ++vv;
	    }
	  
	  //
	  // Swap in any known variables
	  //
	  
	  Helper::swap_in_variables( &l1 , allvars , false , true ); // F = don't allow missing, T = silent                           
	  
	  //
	  // Expand any numerics [][]  -- calls xsigs()
	  //
	  
	  Helper::expand_numerics( &l1 );

	}

      
      //
      // Parse for IF/FI or LOOP/END-LOOP
      //
      
      std::vector<std::string> ctok = Helper::quoted_parse( l1 , "\t " );

      if ( ctok.size() == 0 ) continue;

      const std::string str = Helper::toupper( ctok[0] );
      	
      //
      // Deal with conditionals first: if if_count>1, implies to
      // ignore -- unless we come across an ENDIF/FI
      //
            
      const bool is_if = str == "IF" || str == "IFNOT" ;
      const bool is_fi = str == "ENDIF" || str == "FI" ;     

      if ( is_fi && if_count == 0 )
	Helper::halt( "unbalanced IF/FI conditionals encountered" );

      // get an IF/IFNOT condition
      if ( is_if )
	{
	  if ( ctok.size() !=2 ) 
	    Helper::halt( str + " command expecting a single logical expression: " + l1 );
	  
	  // save condition:

	  // IF ${x}   or   IF #{x}   IF ?{x} 
	  //  this may be a variable/index,
	  //    e.g. which must be set if so, but by here will aleady have been parsed (e.g. to a truth value)

	  if_condition = ctok[1];
	  
	}

      //
      // are we currently inside a spliced-out section?
      //  -- if so, move on unless we have an FI/ENDIF here?
      //
      
      if ( excise )
	{
	  
	  //
	  // Pay attention to any ENDIF/FI statements
	  //
	  
	  if ( is_fi )
	    {
	      if_condition = "";
	      
	      if ( if_count == 0 )
		Helper::halt( "unbalanced IF/FI pairs: " + l1 );

	      --if_count;

	      if ( if_count == 0 )
		excise = false; 
	      
	    }

	  // if we encounter an IF when spliced, still increase the count, but
	  // this has no effect on the splice status
	  if ( is_if )
	    {
	      ++if_count;
	    }
	  
	  // but in all cases, continue to skip these lines due to a prior IF
	  continue;
	  
	}
      else if ( is_if )
	{
	  
	  const bool ifnot = str == "IFNOT";
	  
	  const bool val = Helper::yesno( if_condition );

	  // track depth of IFs
	  ++if_count;
	  
	  if ( ifnot ) // requires F
	    {
	      if ( val ) excise = true;		
	    }
	  else // requires T
	    {
	      if ( ! val ) excise = true;		
	    }
	  
	  continue;
	}
      
      //
      // ignore END/FI for executed blocks
      //
      
      if ( is_fi )
	{
	  if ( if_count == 0 )
	    Helper::halt( "unbalanced IF/FI pairs: " + l1 );
	  
	  --if_count;
	  
	  if ( if_count == 0 )
	    excise = false;
	  
	  continue;
	}

      //
      // Now handle loops
      //

      const bool is_repeat = str == "LOOP";

      const bool is_until  = str == "END-LOOP";
      
      //
      // Are we starting a new loop? 
      //
      
      if ( is_repeat )
	{

	  // track depth of loops
	  ++loop_count;

	  // LOOP index=x vals=1,2,3,4	  
	  param_t param;
          for (int j=1;j<ctok.size();j++)
	    param.parse( ctok[j] );

	  // use separate #{} for loop indexes, as we
	  // are not allowed to assign to them, unlike ${var=value}
	  // pairs; also ${a} != #{a}
	  
	  inner = "#{" + param.requires( "index" ) + "}";
	  
	  // set repeat pointer to start of loop
	  repeat_pointer[ inner ] = i;
	  
	  // note this is the current inner loop (end element)
	  keys.push_back( inner );
	  
	  // at this point, variables and other expansions should have been resolved
	  const std::string expr = param.value( "vals" ) ;

	  // allow inc/exc too?
	  seqs[ inner ] = Helper::parse( Helper::incexc( Helper::xsigs( expr ) ) , "," );
	  
	  if ( seqs[ inner ].size() == 0 )
	    Helper::halt( "empty vals in LOOP for " + inner );

	  // set at first slot
	  slots[ inner ] = 0;

	  // set value
	  vals[ inner ] = seqs[ inner ][ slots[ inner ] ];
	  	  
	}
      

      //
      // Or loop back or ending an existing one?
      //

      if ( is_until )
	{
	  
	  // advance to next slot
          slots[ inner ]++;
	  
	  // all done?
	  const bool closeout = seqs[ inner ].size() == slots[ inner ] ;

	  if ( ! closeout ) // i.e. will be looping back
	    {
	      // update value
	      vals[ inner ] = seqs[ inner ][ slots[ inner ] ];
	      
	      // reset point to read from in input
	      i = repeat_pointer[ inner ];

	    }
	  else // i.e. done w/ the loop
	    {

	      // track loop depth
	      if ( loop_count == 0 )
		Helper::halt( "unbalanced LOOP/END-LOOPs: " + l1 );
	      else
		--loop_count;

	      
	      seqs.erase( seqs.find( inner ) );	      
	      slots.erase( slots.find( inner ) );
	      vals.erase( vals.find( inner ) );
	      repeat_pointer.erase( repeat_pointer.find( inner ) );
	      
	      if ( keys.size() == 1 )
		{
		  inner = "";
		  keys.clear();
		}
	      else
		{
		  keys.pop_back();
		  inner = keys.back();
		}
	    }
	  
	}
      
      
      //
      // Otherwise, just track this line - but swapping in any loop variables
      //
      
      if ( ! ( is_repeat || is_until ) )
	{
	  r.push_back( l1 );	      
	  if ( globals::mirror )
	    logger << "  " << Helper::ltrim( l1 ) << "\n";
	}
      
      // next input line
    }

  if ( globals::mirror )
    logger << "\n ------------------------------------------------------------\n\n";
  
  return r;

}


std::vector<std::string> cmd_t::proc_forloops( const std::vector<std::string> & lines )
{

  std::vector<std::string> r;
  return r;

  
  // // allow multiple nested FOR UNTIL statements
  // int fors = 0 , untils = 0;

  // for (int i=0; i<lines.size(); i++)
  //   {
  //     std::vector<std::string> ctok = Helper::quoted_parse( lines[i] , "\t " );
  //     if ( ctok.size() == 0 ) continue;
  //     const std::string str = Helper::toupper( ctok[0] );
  //     if      ( str == "LOOP" ) ++fors;
  //     else if ( str == "END-LOOP" ) ++untils; 
  //   }
  
  // if ( fors != untils )
  //   Helper::halt( "unbalanced number of LOOP/END-LOOP commands in the script" );

  // // nothing to do? 
  // if ( fors == 0 ) return lines;
  
  // // make expansions
  
  // std::vector<std::string> keys;    // for nested loops, last is the inner
  // std::map<std::string,int> slots;  // for a given loop, which slot we at?
  // std::map<std::string,std::vector<std::string> > seqs; // for a given loop, what are the slots?
  // std::map<std::string,std::string> vals; // swap in this value
  // std::map<std::string,int> repeat_pointer; // line to loop back to 

  // std::string inner = ""; // inner loop
  // const int nl = lines.size();

  // for (int i=0; i<nl; i++)
  //   {

  //     // swap any loop indexes (e.g. as they may feature in a LOOP vals= arg)
      
  //     std::string l1 = lines[i];
      
  //     std::map<std::string,std::string>::const_iterator vv = vals.begin();
  //     while ( vv != vals.end() )
  // 	{
  // 	  l1 = Helper::search_replace( l1 , vv->first , vv->second ) ;
  // 	  ++vv;
  // 	}

  //     // parse
  //     std::vector<std::string> ctok = Helper::quoted_parse( l1 , "\t " );
  //     if ( ctok.size() == 0 ) continue;
  //     const std::string str = Helper::toupper( ctok[0] );
  //     const bool is_repeat = str == "LOOP";
  //     const bool is_until  = str == "END-LOOP";
      
  //     //
  //     // Are we starting a new loop? 
  //     //
      
  //     if ( is_repeat )
  // 	{
	  
  // 	  // LOOP index=x vals=1,2,3,4
	  
  // 	  param_t param;
  //         for (int j=1;j<ctok.size();j++)
  // 	    param.parse( ctok[j] );
	  
  // 	  inner = "#{" + param.requires( "index" ) + "}";
	  
  // 	  // set repeat pointer to start of loop
  // 	  repeat_pointer[ inner ] = i;
	  
  // 	  // note this is the current inner loop (end element)
  // 	  keys.push_back( inner );

  // 	  // allow [1:n] expansion for sequences here
  // 	  const std::string expr = param.value( "vals" ) ;
  // 	  if ( expr.find( "${" ) != std::string::npos )
  // 	    Helper::halt( "currently, LOOP vals cannot contain variables" );	  

  // 	  // allow inc/exc too?
  // 	  seqs[ inner ] = Helper::parse( Helper::incexc( Helper::xsigs( expr ) ) , "," );
	  
  // 	  if ( seqs[ inner ].size() == 0 )
  // 	    Helper::halt( "empty vals in LOOP for " + inner );

  // 	  // set at first slot
  // 	  slots[ inner ] = 0;

  // 	  // set value
  // 	  vals[ inner ] = seqs[ inner ][ slots[ inner ] ];
	  	  
  // 	}
      
  //     //
  //     // Or loop back or ending an existing one?
  //     //

  //     if ( is_until )
  // 	{
	  
  // 	  // advance to next slot
  //         slots[ inner ]++;
	  
  // 	  // all done?
  // 	  const bool closeout = seqs[ inner ].size() == slots[ inner ] ;

  // 	  if ( ! closeout )
  // 	    {
  // 	      // update value
  // 	      vals[ inner ] = seqs[ inner ][ slots[ inner ] ];
	      
  // 	      // reset point to read from in input
  // 	      i = repeat_pointer[ inner ];

  // 	    }
  // 	  else
  // 	    {
  // 	      seqs.erase( seqs.find( inner ) );	      
  // 	      slots.erase( slots.find( inner ) );
  // 	      vals.erase( vals.find( inner ) );
  // 	      repeat_pointer.erase( repeat_pointer.find( inner ) );
	      
  // 	      if ( keys.size() == 1 )
  // 		{
  // 		  inner = "";
  // 		  keys.clear();
  // 		}
  // 	      else
  // 		{
  // 		  keys.pop_back();
  // 		  inner = keys.back();
  // 		}
  // 	    }
	  
  // 	}

      
  //     //
  //     // Otherwise, just track this line - but have swapped in any loop variables
  //     //
      
  //     if ( ! ( is_repeat || is_until ) )
  // 	{
  // 	  r.push_back( l1 );	      
  // 	}
  //     // next input line
  //   }

  
  // return r;
}

void cmd_t::dump()
{
  for (int i=0; i<cmds.size(); i++)
    {      
      std::cout << cmds[i] << "\n";
      param_t param = params[i];
      std::cout << param.dump() << "\n\n";
    }
}

std::string cmd_t::concat( const std::vector<std::string> & l )
{
  std::string s;
  for (int i=0; i<l.size(); i++)
    s += l[i] + "\n";
  return s;
}
