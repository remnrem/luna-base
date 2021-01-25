
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
#include "luna.h"

extern logger_t logger;
extern writer_t writer;


//
// param_t 
//

void param_t::add( const std::string & option , const std::string & value ) 
{
  
  if ( opt.find( option ) != opt.end() ) 
    Helper::halt( option + " parameter specified twice, only one value would be retained" );
  
    opt[ option ] = value; 
}  


void param_t::add_hidden( const std::string & option , const std::string & value ) 
{
  add( option , value );
  hidden.insert( option );
}

int param_t::size() const 
{ 
  // handle hidden things...
  return opt.size() - hidden.size();
}

void param_t::parse( const std::string & s )
{
  std::vector<std::string> tok = Helper::quoted_parse( s , "=" );
  if ( tok.size() == 2 )     add( tok[0] , tok[1] );
  else if ( tok.size() == 1 ) add( tok[0] , "T" );
  else // ignore subsequent '=' signs in 'value'  (i.e. key=value=2  is 'okay', means "value=2" is set to 'key')
    {
      std::string v = tok[1];
      for (int i=2;i<tok.size();i++) v += "=" + tok[i];
      add( tok[0] , v );
    }
}

void param_t::update( const std::string & id , const std::string & wc )
{

  // replace all instances of 'globals::indiv_wildcard' with 'id'
  // for all values
  
  std::map<std::string,std::string>::iterator ii = opt.begin();
  while ( ii != opt.end() ) 
    {

      std::string v = ii->second;
      bool changed = false;

      // 1. replace indiv wildcard (e.g. ^) with this person's ID
      //    nb.  this will also happen via ${id} which is a special, automatic individual-level
      //         variable
      
      while ( v.find( wc ) != std::string::npos )
	{
	  int p = v.find( wc );
	  v = v.substr( 0 , p ) + id + v.substr(p+1);
	  changed = true;
	}
      
      // 2. for any @{includes}, insert contents of file (comma-delimited)
      
      if ( Helper::swap_in_includes( &v ) )
	changed = true;

      //
      // needs updating?
      //

      if ( changed || v != ii->second ) 
	ii->second = v;
      
      ++ii;
    }
  
  
}

void param_t::clear() 
{ 
  opt.clear(); 
  hidden.clear(); 
} 

bool param_t::has(const std::string & s ) const 
{
  return opt.find(s) != opt.end(); 
} 

std::string param_t::value( const std::string & s ) const 
{ 
  if ( has( s ) )
    return opt.find( s )->second;
  else
    return "";
}

bool param_t::single() const 
{ 
  return size() == 1; 
}

std::string param_t::single_value() const 
{ 
  if ( ! single() ) Helper::halt( "no single value" ); 
  
  std::map<std::string,std::string>::const_iterator ii = opt.begin();
      
  while ( ii != opt.end() ) 
    {
      if ( hidden.find( ii->first ) == hidden.end() ) return ii->first;
      ++ii;
    }
  return ""; // should not happen
}

std::string param_t::requires( const std::string & s ) const
{
  if ( ! has(s) ) Helper::halt( "command requires parameter " + s );
  return value(s);
}

int param_t::requires_int( const std::string & s ) const
{
  if ( ! has(s) ) Helper::halt( "command requires parameter " + s );
  int r;
  if ( ! Helper::str2int( value(s) , &r ) ) 
    Helper::halt( "command requires parameter " + s + " to have an integer value" );
  return r;
}

double param_t::requires_dbl( const std::string & s ) const
{
  if ( ! has(s) ) Helper::halt( "command requires parameter " + s );
  double r;
  if ( ! Helper::str2dbl( value(s) , &r ) ) 
    Helper::halt( "command requires parameter " + s + " to have a numeric value" );
  return r;
}

std::string param_t::dump( const std::string & indent , const std::string & delim ) const
{
  std::map<std::string,std::string>::const_iterator ii = opt.begin();
  int sz = opt.size();
  int cnt = 1;
  std::stringstream ss;
  while ( ii != opt.end() ) 
    {
      if ( cnt == sz )
	ss << indent << ii->first << "=" << ii->second; 
      else
	ss << indent << ii->first << "=" << ii->second << delim; 
      ++cnt;
      ++ii;
    }
  return ss.str();
}

std::set<std::string> param_t::strset( const std::string & k , const std::string delim ) const
{
  std::set<std::string> s;
  if ( ! has(k) ) return s;
  std::vector<std::string> tok = Helper::quoted_parse( value(k) , delim );
  for (int i=0;i<tok.size();i++) s.insert( Helper::unquote( tok[i]) );
  return s;
}

std::vector<std::string> param_t::strvector( const std::string & k , const std::string delim ) const
{
  std::vector<std::string> s;
  if ( ! has(k) ) return s;
  std::vector<std::string> tok = Helper::quoted_parse( value(k) , delim );
  for (int i=0;i<tok.size();i++) s.push_back( Helper::unquote( tok[i]) );
  return s;
}

std::vector<double> param_t::dblvector( const std::string & k , const std::string delim ) const
{
  std::vector<double> s;
  if ( ! has(k) ) return s;
  std::vector<std::string> tok = Helper::quoted_parse( value(k) , delim );
  for (int i=0;i<tok.size();i++) 
    {
      std::string str = Helper::unquote( tok[i]);
      double d = 0;
      if ( ! Helper::str2dbl( str , &d ) ) Helper::halt( "Option " + k + " requires a double value(s)" );
      s.push_back(d); 
    }
  return s;
}

std::vector<int> param_t::intvector( const std::string & k , const std::string delim ) const
{
  std::vector<int> s;
  if ( ! has(k) ) return s;
  std::vector<std::string> tok = Helper::quoted_parse( value(k) , delim );
  for (int i=0;i<tok.size();i++) 
    {
      std::string str = Helper::unquote( tok[i]);
      int d = 0;
      if ( ! Helper::str2int( str , &d ) ) Helper::halt( "Option " + k + " requires an integer value(s)" );
      s.push_back(d);
    }
  return s;
}


std::set<std::string> param_t::keys() const
{
  std::set<std::string> s;
  std::map<std::string,std::string>::const_iterator ii = opt.begin();
  while ( ii != opt.end() )
    {
      s.insert( ii->first );
      ++ii;
    }
  return s;
}



//
// cmd_t
//


cmd_t::cmd_t() 
{
  register_specials();
  reset();
  error = ! read();
}

cmd_t::cmd_t( const std::string & str ) 
{
  register_specials();
  reset();
  error = ! read( &str , true ); 
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
  
  // however: two rules
  // 1. many-to-one mapping means the same label cannot have multiple primary aliases
  // 2. following, and on principle of no transitive properties, alias cannot have alias

  // X|Y|Z
  // X|A|B
  
  // W|A  bad, A already mapped
  // V|X  bad, X already mapped
  // i.e. things can only occur once in the RHS, or multiple times in the LHS
  
  
  // format canonical|alias1|alias2 , etc.
  std::vector<std::string> tok = Helper::quoted_parse( s , "|" );    
  if ( tok.size() < 2 ) Helper::halt( "bad format for signal alias:  canonical|alias 1|alias 2\n" + s );
  const std::string primary = Helper::unquote( tok[0] );
  for (int j=1;j<tok.size();j++) 
    {

      // impose rules
      const std::string mapped = Helper::unquote( tok[j] ) ;
      
      if ( primary_alias.find( mapped ) != primary_alias.end() )
	Helper::halt( mapped + " specified as both primary alias and mapped term" );

      if ( label_aliases.find( mapped ) != label_aliases.end() )
	if ( primary != label_aliases[ mapped ] )  
	  Helper::halt( mapped + " specified twice in alias file w/ different primary aliases" );

      // otherwise, set 
      label_aliases[ mapped ] = primary;

      primary_alias[ primary ].push_back( mapped );
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

void cmd_t::replace_wildcards( const std::string & id )
{
  
  //
  // get copy of original script;  comments, newlines vesrus '&' will already have
  // been taken care of
  //
  
  std::string iline = line;


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
  //
  // remove any conditional blocks, where variable should be var=1 (in) or var=0 (out)
  // or a value for add=variable,var2 for var=1 var=2
  //

  // [[var1
  //   block
  // ]]var1
  
  Helper::process_block_conditionals( &iline , allvars );

  
  //
  // swap in any variables (and allows for them being defined on-the-fly)
  //

  Helper::swap_in_variables( &iline , &allvars );

  //
  // expand [var][1:10] sequences
  //

  Helper::expand_numerics( &iline );


  //  std::cerr << "final [" << iline << "]\n";
  
  //
  // Parse into commands/options
  //

  std::vector<std::string> tok = Helper::quoted_parse( iline , "\n" );
  
  // command(s): do this just for printing; real parsing will include variables, etc

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
  
  //params = original_params;
  for (int p = 0 ; p < params.size(); p++ )
    params[p].update( id , globals::indiv_wildcard );

   
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
	  
	  // is this a continuation line?
	  bool continuation = s[0] == ' ' || s[0] == '\t';

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
    allinput.str( cmd_t::cmdline_cmds );
      
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
  // Note... used to do processing of command file here (global)
  //  as we now use ivars as well as vars, hold on this till later,
  //  in update()
  //


  //
  // Initial processing 
  //
  
  std::vector<std::string> tok = Helper::quoted_parse( line , "\n" );
  
  if ( tok.size() == 0 ) 
    {
      quit(true);      
      return false;
    }


  
  // command(s): do this just for printing; real parsing will include variables, etc
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


  return true;
}



//
// Evaluate commands
//

bool cmd_t::eval( edf_t & edf ) 
{

  //
  // Loop over each command
  //
  
  for ( int c = 0 ; c < num_cmds() ; c++ )
    {	        

      // was a problem flag raised when loading the EDF?
      
      if ( globals::problem ) return false;
      
      //
      // If this particular command did not explicitly specify
      // signals, then add all
      //
      
      if ( ! param(c).has( "sig" ) )
	param(c).add_hidden( "sig" , signal_string() );
      

      //
      // Print command
      //

      logger << " ..................................................................\n"
	     << " CMD #" << c+1 << ": " << cmd(c) << "\n";
      logger << "   options: " << param(c).dump( "" , " " ) << "\n";
      
      
      writer.cmd( cmd(c) , c+1 , param(c).dump( "" , " " ) );

      // use strata to keep track of tables by commands, with leading underscore to denote
      // special status of this factor

      writer.level( cmd(c) , "_" + cmd(c) );
      
      
      //
      // Now process the command
      //
      
      if      ( is( c, "WRITE" ) )        proc_write( edf, param(c) );
      else if ( is( c, "SUMMARY" ) )      proc_summaries( edf , param(c) );
      else if ( is( c, "HEADERS" ) )      proc_headers( edf , param(c) );
      else if ( is( c, "ALIASES" ) )      proc_aliases( edf , param(c) );

      else if ( is( c, "DESC" ) )         proc_desc( edf , param(c) );
      else if ( is( c, "TYPES" ) )        proc_show_channel_map();
      else if ( is( c, "VARS" ) )         proc_dump_vars( edf , param(c) );
      else if ( is( c, "STATS" ) )        proc_stats( edf , param(c) );
      
      else if ( is( c, "REFERENCE" ) )    proc_reference( edf , param(c) );
      else if ( is( c, "DEREFERENCE" ) )  proc_dereference( edf , param(c) );
      else if ( is( c, "ADJUST" ) )       proc_adjust( edf , param(c) ); 

      else if ( is( c, "FLIP" ) )         proc_flip( edf , param(c) );
      else if ( is( c, "CANONICAL" ) )    proc_canonical( edf , param(c) );
      else if ( is( c, "uV" ) )           proc_scale( edf , param(c) , "uV" ); 
      else if ( is( c, "mV" ) )           proc_scale( edf , param(c) , "mV" );
      else if ( is( c, "MINMAX" ) )       proc_minmax( edf , param(c) );
      
      else if ( is( c, "RECORD-SIZE" ) )  proc_rerecord( edf , param(c) );
      
      else if ( is( c, "TIME-TRACK" ) )   proc_timetrack( edf, param(c) );
      //      else if ( is( c, "CONTIN" ) )       proc_continuous( edf, param(c) ); // not allowed....

      else if ( is( c, "STAGE" ) )        proc_sleep_stage( edf , param(c) , false );
      else if ( is( c, "HYPNO" ) )        proc_sleep_stage( edf , param(c) , true );
      
      else if ( is( c, "TSLIB" ) )        pdc_t::construct_tslib( edf , param(c) );
      else if ( is( c, "SSS" ) )          pdc_t::simple_sleep_scorer( edf , param(c) );
      else if ( is( c, "EXE" ) )          pdc_t::similarity_matrix( edf , param(c) );
      //      else if ( is( c, "CHCHK" ) )        pdc_t::channel_checker( edf , param(c) );

      //      else if ( is( c, "LW" ) )           lw_prep_t( edf , param(c) );
      
      else if ( is( c, "DUMP" ) )         proc_dump( edf, param(c) );      
      else if ( is( c, "DUMP-RECORDS" ) ) proc_record_dump( edf , param(c) );
      else if ( is( c, "RECS" ) )         proc_record_table( edf , param(c) );
      else if ( is( c, "SEGMENTS" ) )     proc_dump_segs( edf , param(c) );
      else if ( is( c, "DUMP-EPOCHS" ) )  proc_epoch_dump( edf, param(c) ); // REDUNDANT; use ANNOTS epoch instead

      else if ( is( c, "ANNOTS" ) )       proc_list_all_annots( edf, param(c) );
      else if ( is( c, "WRITE-ANNOTS" ) ) proc_write_annots( edf, param(c) );
      else if ( is( c, "A2S" ) )          proc_annot2signal( edf, param(c) );
      else if ( is( c, "SPANNING" ) ) proc_list_spanning_annots( edf, param(c) );
      //else if ( is( c, "COUNT-ANNOTS" ) ) proc_list_annots( edf , param(c) ); // REDUNDANT; use ANNOTS epoch instead

      else if ( is( c, "MATRIX" ) )       proc_epoch_matrix( edf , param(c) );
      else if ( is( c, "RESTRUCTURE" ) || is( c, "RE" ) )  proc_restructure( edf , param(c) );
      else if ( is( c, "SIGNALS" ) )      proc_drop_signals( edf , param(c) );
      else if ( is( c, "COPY" ) )         proc_copy_signal( edf , param(c) );
      else if ( is( c, "ORDER" ) )        proc_order_signals( edf , param(c) );

      else if ( is( c, "RMS" ) || is( c, "SIGSTATS" ) ) proc_rms( edf, param(c) );
      else if ( is( c, "MSE" ) )          proc_mse( edf, param(c) );
      else if ( is( c, "LZW" ) )          proc_lzw( edf, param(c) );
      else if ( is( c, "ZR" ) )           proc_zratio( edf , param(c) );
      else if ( is( c, "ANON" ) )         proc_anon( edf , param(c) );
      else if ( is( c ,"EPOCH" ) )        proc_epoch( edf, param(c) );
      else if ( is( c ,"SLICE" ) )        proc_slice( edf , param(c) , 1 );

      else if ( is( c , "SUDS" ) )        proc_suds( edf , param(c) );
      else if ( is( c , "MAKE-SUDS" ) )   proc_make_suds( edf , param(c) );
      else if ( is( c , "SOAP-SUDS" ) )   proc_self_suds( edf , param(c) );

      else if ( is( c, "EVAL" ) )         proc_eval( edf, param(c) );
      else if ( is( c, "MASK" ) )         proc_mask( edf, param(c) );

      else if ( is( c, "FILE-MASK" ) )    proc_file_mask( edf , param(c) ); // not supported/implemented
      else if ( is( c, "DUMP-MASK" ) )    proc_dump_mask( edf, param(c) );
      else if ( is( c, "CHEP" ) )         timeline_t::proc_chep( edf, param(c) );
      else if ( is( c, "CHEP-MASK" ) )    proc_chep_mask( edf, param(c) );
      
      else if ( is( c, "EPOCH-ANNOT" ) )  proc_file_annot( edf , param(c) );
      else if ( is( c, "EPOCH-MASK" ) )   proc_epoch_mask( edf, param(c) );
      
      
      else if ( is( c, "FILTER" ) )       proc_filter( edf, param(c) );      
      else if ( is( c, "FILTER-DESIGN" )) proc_filter_design( edf, param(c) );
      else if ( is( c, "CWT-DESIGN" ) )   proc_cwt_design( edf , param(c) );
      else if ( is( c, "CWT" ) )          proc_cwt( edf , param(c) );
      else if ( is( c, "HILBERT" ) )      proc_hilbert( edf , param(c) );

      else if ( is( c, "TV" ) )           proc_tv_denoise( edf , param(c) );
      
      else if ( is( c, "COVAR" ) )        proc_covar( edf, param(c) );
      else if ( is( c, "PSD" ) )          proc_psd( edf, param(c) );	  
      else if ( is( c, "MTM" ) )          proc_mtm( edf, param(c) );
      else if ( is( c, "1FNORM" ) )       proc_1overf_norm( edf, param(c) );

      else if ( is( c, "PSC" ) )          proc_psc( edf , param(c) );
      
      else if ( is( c, "MS" ) )           proc_microstates( edf , param(c) );

      else if ( is( c, "TLOCK" ) )        proc_tlock( edf , param(c) );
      else if ( is( c, "PEAKS" ) )        proc_peaks( edf , param(c) );

      else if ( is( c, "SEDF" ) )         proc_sedf( edf , param(c) );

      else if ( is( c, "FIP" ) )          proc_fiplot( edf , param(c) );
      
      else if ( is( c, "COH" ) )          proc_coh( edf , param(c) );
      else if ( is( c, "CC" ) )           proc_conncoupl( edf , param(c) );
      else if ( is( c, "CORREL" ) )       proc_correl( edf , param(c) );
      else if ( is( c, "PSI" ) )          proc_psi( edf , param(c) );
      else if ( is( c, "ACF" ) )          proc_acf( edf , param(c) );
      else if ( is( c, "ED" ) )           proc_elec_distance( edf , param(c) );
      else if ( is( c, "ICA" ) )          proc_ica( edf, param(c) );
      else if ( is( c, "CLOCS" ) )        proc_attach_clocs( edf , param(c) );
      else if ( is( c, "L1OUT" ) )        proc_leave_one_out( edf , param(c) );
      else if ( is( c, "INTERPOLATE" ) )  proc_chep_based_interpolation( edf, param(c) );
      else if ( is( c, "SL" ) )           proc_surface_laplacian( edf , param(c) );
      else if ( is( c, "EMD" ) )          proc_emd( edf , param(c) );
      
      else if ( is( c, "MI" ) )           proc_mi( edf, param(c) );
      else if ( is( c, "HR" ) )           proc_bpm( edf , param(c) );
      else if ( is( c, "SUPPRESS-ECG" ) ) proc_ecgsuppression( edf , param(c) );
      else if ( is( c, "PAC" ) )          proc_pac( edf , param(c) );
      else if ( is( c, "CFC" ) )          proc_cfc( edf , param(c) );
      else if ( is( c, "TAG" ) )          proc_tag( param(c) );
      else if ( is( c, "RESAMPLE" ) )     proc_resample( edf, param(c) );
      
      else if ( is( c, "SPINDLES" ) )     proc_spindles( edf, param(c) );	  
      else if ( is( c, "SO" ) )           proc_slowwaves( edf, param(c) );
      else if ( is( c, "COUPL" ) )        proc_coupling( edf , param(c) );
      
      else if ( is( c, "POL" ) )          proc_polarity( edf, param(c) );	  
      else if ( is( c, "REMS" ) )         proc_rems( edf, param(c) );
      
      else if ( is( c, "ARTIFACTS" ) )    proc_artifacts( edf, param(c) );

      else if ( is( c, "CACHE" ) )        proc_dump_cache( edf , param(c) );
      
      else if ( is( c, "SPIKE" ) )        proc_spike( edf , param(c) );
      else if ( is( c, "SHIFT" ) )        proc_shift( edf , param(c) );
      
      else 
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
// ----------------------------------------------------------------------------------------
//


// HEADERS : summarize EDF files

void proc_headers( edf_t & edf , param_t & param )
{
  // optionally add a SIGNALS col that has a comma-delimited list of all signals
  edf.terse_summary( param );
}

// ALIASES : report aliasing of channels and annotations

void proc_aliases( edf_t & edf , param_t & param )
{
  edf.report_aliases();
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



// SOAP-SUDS : single observation accuracies and probabilities;
//   i.e. use SUDS on self, to evaluate staging/signal quality

void proc_self_suds( edf_t & edf , param_t & param  )
{  
  suds_t::set_options( param );  
  suds_indiv_t self;
  self.evaluate( edf , param );  
}

// MAKE-SUDS : populate folder 'db' with trainers

void proc_make_suds( edf_t & edf , param_t & param  )
{  
  suds_t::set_options( param );  
  suds_indiv_t trainer;
  trainer.add_trainer( edf , param );  
}


// SUDS : stageing

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
  
  // this is only done once per session, i.e. even if multiple targets
  // are scored

  // bank() and wbank() can share the same individuals (they will only
  // be loaded once) load wbank() first, as that also involves loading
  // the PSD (i.e. raw features). If the an individual is only in the
  // bank, these are not needed/loaded

  //
  // Weight trainers
  //

  // by default, use self-trainer as the weight trainer; if param
  // 'wdb' is given explicitly, then ALL indivs in that database will be
  // used to retrain the trainer weights
  
  if ( param.has( "wdb" ) ) 
    suds.attach_db( param.value( "wdb" ) , true );
  else // else, use the same training panel (but also loading raw features)
    suds.attach_db( param.value( "db" ) , true );

 
  //
  // Trainers (if not already loaded by wbank; in this case,
  // attach_db() will just skip that person
  //

  // attaach_db() F -> do not load PSD (not needed)
  suds.attach_db( param.requires( "db" ) , false );

  
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


// ARTIFACTS : artifact rejection using Buckelmueller et al. 

void proc_artifacts( edf_t & edf , param_t & param )	  
{
  std::string signal = param.requires( "sig" );
  annot_t * a = buckelmuller_artifact_detection( edf , param , signal );  
}


// LEGACY-FILTER : band-pass filter, band-pass filter only

// void proc_filter_legacy( edf_t & edf , param_t & param )	  
// {
//   //band_pass_filter( edf , param );
// }


// FILTER : general FIR

void proc_filter( edf_t & edf , param_t & param )	  
{
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


// RESAMPLE : band-pass filter

void proc_resample( edf_t & edf , param_t & param ) 
{
  dsptools::resample_channel( edf, param );
}


// MS: microstate analysis
void proc_microstates( edf_t & edf , param_t & param )
{
  dsptools::microstates( edf , param );
}


// TLOCK
void proc_tlock( edf_t & edf  , param_t & param )
{
  // get mean time-locked value of one signal against a set of annotations (time-points)
  dsptools::tlock( edf , param );
}

// PEAKS
void proc_peaks( edf_t & edf , param_t & param )
{
  dsptools::peaks( edf , param );
}

// SEDF : make a summarize EDF 
void proc_sedf( edf_t & edf , param_t & param )
{
  sedf_t sedf( edf , param );  
}



// PSC : either build PSC (from multiple results) or fit to an EDF
void proc_psc( edf_t & edf , param_t & param )
{

  // note:  PSC contruct() is called w/out an EDF from the command line
  //  luna --psc <args>

  psc_t psc;
  
  // if already populated, this returns so we can call multiple times
  // expects 'proj=' variable to point to output of 'proj=' from PSC 
  psc.attach( param );
  
  // project attached solution
  psc.project( edf , param );

}

// PSD : calculate PSD 

void proc_psd( edf_t & edf , param_t & param )	  
{  
  std::string signal = param.requires( "sig" );
  annot_t * power = spectral_power( edf , signal , param );  
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

// FI-plot : frequency/interval plot

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
    logger << " setting analysis tag to [" << globals::current_tag << "]\n";

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
      
      edf.header.patient_id = edf.id;
      
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


// SPINDLES : spindle detection using CWT or bandpass/RMS

void proc_spindles( edf_t & edf , param_t & param )
{	

  // default to wavelet analysis
  std::string method = param.has( "method" ) ? param.value( "method" ) : "wavelet" ; 
  
  annot_t * a = NULL;

  if      ( method == "bandpass" ) a = spindle_bandpass( edf , param );
  else if ( method == "wavelet" ) a = spindle_wavelet( edf , param );
  else Helper::halt( "SPINDLE method not recognized; should be 'bandpass' or 'wavelet'" );

}

// COUPL : spindle/SO couplig

void proc_coupling( edf_t & edf , param_t & param )
{
  // requires cached SPINDLES and SO results
  spindle_so_coupling( edf , param );
}


// POL : polarity check for EEG N2/N3 

void proc_polarity( edf_t & edf , param_t & param )
{	
  dsptools::polarity( edf , param );
}

// REMS : detect REMS via simple heuristic

void proc_rems( edf_t & edf , param_t & param )
{
  dsptools::rems( edf , param );
}

// SW || SLOW-WAVES : detect slow waves, do time-locked FFT on rest of signal

void proc_slowwaves( edf_t & edf , param_t & param )
{	

  // find slow-waves
  slow_waves_t sw( edf , param );
  
}


// WRITE : write a new EDF or EDFZ to disk
// optionally, writes annotation files,
// (although user responsible for not altering time structure of EDF
// i.e. we do not change the time encoding of the annotation file, which
// are always anchored to the original)

void proc_write( edf_t & edf , param_t & param )
{
  
  // write a .edfz and .edfz.idx
  bool edfz = param.has( "edfz" );
    
  // add 'tag' to new EDF
  std::string filename = edf.filename;

  if ( Helper::file_extension( filename, "edf" ) || 
       Helper::file_extension( filename, "EDF" ) ) 
    filename = filename.substr(0 , filename.size() - 4 );
  
  if ( Helper::file_extension( filename, "edfz" ) || 
       Helper::file_extension( filename, "EDFZ" ) ) 
    filename = filename.substr(0 , filename.size() - 5 );

  filename += "-" + param.requires( "edf-tag" ) + ".edf";
  if ( edfz ) filename += "z";
  
  //
  // optionally, allow directory change
  //

  if ( param.has( "edf-dir" ) )
    {
      const std::string outdir = param.value("edf-dir");
      
      if ( outdir[ outdir.size() - 1 ] != globals::folder_delimiter ) 
	Helper::halt("edf-dir value must end in '" + std::string(1,globals::folder_delimiter) + "' to specify a folder" );

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
      logger << " appending " << filename << " to sample-list " << file << ( append_annots ? " (with annotations)" : " (dropping any annotations)" ) << "\n";
      
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

  bool write_as_edf = param.has( "force-edf" );
  
  if ( write_as_edf ) edf.set_edf();
  
  //
  // Save data (write_as_edf flag forces starttime to 00.00.00)
  //

  bool saved = edf.write( filename , edfz , write_as_edf );

  if ( saved ) 
    logger << " saved new EDF" << ( edf.header.edfplus ? "+" : "" ) << ", " << filename << "\n";

}


// EPOCH : set epoch 

void proc_epoch( edf_t & edf , param_t & param )
{

  // unepoch?
  if ( param.has( "clear" ) )
    {
      logger << "  clearing all epochs: signals are now unepoched\n";
      edf.timeline.unepoch();
      return;
    }
  
  double dur = 0 , inc = 0;
  
  // default = 30 seconds, non-overlapping
  if ( ! ( param.has( "len" ) || param.has("dur") || param.has( "epoch" ) ) )
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


  // if already EPOCH'ed for a different record size, or increment,
  // then we should first remove all epochs; this will remove the
  // EPOCH-to-record mapping too
  
  if ( edf.timeline.epoched() 
       && ( ( ! Helper::similar( edf.timeline.epoch_length() , dur ) )
	    || ( ! Helper::similar( edf.timeline.epoch_inc() , inc ) ) ) )
    {
      logger << " epoch definitions have changed: original epoch mappings will be lost\n";
      edf.timeline.unepoch();
    }

  //
  // basic log info
  //

  
  int ne = edf.timeline.set_epoch( dur , inc );  


  //
  // minimal output to stdout
  //

  if ( param.has( "min" ) )
    {
      std::cout << ne << "\n";
      return;
    }
  
  logger << " set epochs, length " << dur << " (step " << inc << "), " << ne << " epochs\n";

  writer.value( "NE" , ne );
  writer.value( "DUR" , dur );
  writer.value( "INC" , inc );
  

  //
  // write more verbose information to db
  //
  
  if ( param.has( "verbose" ) )
    {

      // track clock time

      clocktime_t starttime( edf.header.starttime );
      
      bool hms = starttime.valid;
      
      
      edf.timeline.first_epoch();
      
      while ( 1 ) 
	{
	  
	  int epoch = edf.timeline.next_epoch();      
	  
	  if ( epoch == -1 ) break;
	  
	  interval_t interval = edf.timeline.epoch( epoch );
	  
	  //      std::cout << "epoch " << epoch << "\t" << interval.as_string() << "\n";
	  
	  // original encoding (i.e. to allows epochs to be connected after the fact
	  writer.epoch( edf.timeline.display_epoch( epoch ) );
	  
	  // if present, original encoding
	  writer.value( "E1" , epoch+1 );
	  
	  writer.value( "INTERVAL" , interval.as_string() );      
	  writer.value( "START"    , interval.start_sec() );
	  writer.value( "MID"      , interval.mid_sec() );
	  writer.value( "STOP"     , interval.stop_sec() );
	  
	  // original time-points
	  
	  if ( hms )
	    {
	      const double sec0 = interval.start * globals::tp_duration;
	      clocktime_t present = starttime;
	      present.advance( sec0 / 3600.0 );
	      std::string clocktime = present.as_string();
	      writer.value( "HMS" , clocktime );
	    }
	
	  
	}		  
      
  
      writer.unepoch();
      
    }


  //
  // any constraints on the min. number of epochs required? 
  //

  if ( param.has("require") )
    {
      int r = param.requires_int( "require" );
      if ( ne < r ) 
	{	  
	  logger << " ** warning for "  
		 << edf.id << " when setting EPOCH: "
		 << "required=" << r << "\t"
		 << "but observed=" << ne << "\n";
	  globals::problem = true;
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
  
  // COMMAND NOT SUPPORTED

  std::set<std::string> vars;
  std::string onelabel;
  
  if ( param.has( "if" ) ) 
    {    
      if ( param.has( "ifnot" ) ) Helper::halt( "both if & ifnot specified" );
      vars = param.strset( "if" );
      onelabel = param.value("if");
      logger << " masking epochs that match " << onelabel << "\n";
    }
  else if ( param.has( "ifnot" ) ) 
    {
      vars = param.strset( "ifnot" );
      onelabel = param.value("ifnot");
      logger << " masking epochs that do not match " << onelabel << "\n";
    }
  else
    Helper::halt( "no if/ifnot specified" );
 
  edf.timeline.apply_simple_epoch_mask( vars , onelabel , param.has("if") );  

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


// DUMP-MASK : output the current mask as an .annot file

void proc_dump_mask( edf_t & edf , param_t & param )
{
  if ( ! param.has("tag") )
    {
      edf.timeline.dumpmask();
      return;
    }

  // otherwise, create an ANNOT file from the MASK, i.e. for viewing
  std::string tag = param.requires( "tag" );
  std::string path = param.has( "path" ) ? param.value("path") : ".";
  bool no_id = ! param.has( "no-id" );
  edf.timeline.mask2annot( path, tag , no_id ); 
}

// CHEP : dump, or convert from CHEP->MASK
// in edf/chep.cpp


// COUNT-ANNOTS : show all annotations for the EDF // REDUNDANT

void proc_list_annots( edf_t & edf , param_t & param )
{
  summarize_annotations( edf , param );
}



// WRITE-ANNOTS : write all annots to disk

void proc_write_annots( edf_t & edf , param_t & param )
{
  edf.timeline.annotations.write( param.requires( "file" ) , param );
}

// A2S : make signbal from ANNOTS

void proc_annot2signal( edf_t & edf , param_t & param )
{
  edf.timeline.annot2signal( param );
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
  edf.add_continuous_time_track();
}

// RESTRUCTURE : flush masked records

void proc_restructure( edf_t & edf , param_t & param )
{
  // just drop MASK'ed records, then reset mask
  edf.restructure( );
}


// DUMP-RECORDS : show all records (i.e. raw data)

void proc_record_dump( edf_t & edf , param_t & param )
{
  edf.add_continuous_time_track();  
  edf.record_dumper( param );
}


// SEGMENTS : show all contiguous segments

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
// HYPNO : verbose report on sleep STAGES     (verbose = F)

void proc_sleep_stage( edf_t & edf , param_t & param , bool verbose )
{
  
  std::string wake   = param.has( "wake" )  ? param.value("wake")  : "" ; 
  std::string nrem1  = param.has( "N1" ) ? param.value("N1") : "" ; 
  std::string nrem2  = param.has( "N2" ) ? param.value("N2") : "" ; 
  std::string nrem3  = param.has( "N3" ) ? param.value("N3") : "" ; 
  std::string nrem4  = param.has( "N4" ) ? param.value("N4") : "" ; 
  std::string rem    = param.has( "REM" )  ? param.value("REM")  : "" ; 
  std::string misc   = param.has( "?" )  ? param.value("?")  : "" ; 

  // either read these from a file, or display
  
  if ( param.has( "file" ) )
    {
      std::vector<std::string> ss = Helper::file2strvector( param.value( "file" ) );
      edf.timeline.hypnogram.construct( &edf.timeline , param , verbose , ss );
    }
  else
    {      
      edf.timeline.annotations.make_sleep_stage( wake , nrem1 , nrem2 , nrem3 , nrem4 , rem , misc );
      bool okay = edf.timeline.hypnogram.construct( &edf.timeline , param , verbose ); 
      if ( ! okay ) return; // i.e. if no valid annotations found
    }

  // and output...
  edf.timeline.hypnogram.output( verbose );

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

void proc_attach_clocs( edf_t & edf , param_t & param )
{
  
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

// ICA : fastICA on sample by channel matrix (whole trace)

void proc_ica( edf_t & edf , param_t & param )
{
  dsptools::ica_wrapper( edf , param );
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


// CACHE : internal command to dump cache contents (debugging)

void proc_dump_cache( edf_t & edf , param_t & param )
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
      std::cout << cache->print();
    }
  else if ( str_cache )
    {
      cache_t<std::string> * cache = edf.timeline.cache.find_str( cname );
      if ( cache == NULL ) Helper::halt( "could not find str-cache " + cname );
      std::cout << cache->print();
    }  
  else if ( num_cache )
    {
      cache_t<double> * cache = edf.timeline.cache.find_num( cname );
      if ( cache == NULL ) Helper::halt( "could not find num-cache " + cname );
      std::cout << cache->print();
    }
  else if ( tp_cache )
    {
      cache_t<uint64_t> * cache = edf.timeline.cache.find_tp( cname );
      if ( cache == NULL ) Helper::halt( "could not find tp-cache " + cname );
      std::cout << cache->print();
    }
  
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

// CFC : generic cross-frequency coupling methods (other than PAC as above)

void proc_cfc( edf_t & edf , param_t & param )
{
  dsptools::cfc( edf , param );
}

// SUPPRESS-ECG : ECG supression

void proc_ecgsuppression( edf_t & edf , param_t & param )
{
  dsptools::ecgsuppression( edf , param );
}

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


// COPY : mirror a signal
void proc_copy_signal( edf_t & edf , param_t & param )
{
  
  signal_list_t originals = edf.header.signal_list( param.requires( "sig" ) );

  std::string tag = param.requires( "tag" );
  
  for (int s=0;s<originals.size();s++)
    {

      if ( edf.header.is_data_channel( originals(s) ) )
	{

	  std::string new_label = originals.label( s ) + "_" + tag; 
	  
	  if ( ! edf.header.has_signal( new_label ) )
	    {
	      logger << " copying " << originals.label(s) << " to " << new_label << "\n";
	      edf.copy_signal( originals.label(s) , new_label );
	    }
	}
    }
}

// DROP : drop a signal

void proc_drop_signals( edf_t & edf , param_t & param )
{
  
  std::set<std::string> keeps, drops;
  if ( param.has( "keep" ) ) keeps = param.strset( "keep" );

  if ( param.has( "keep" ) && param.has( "req" ) ) 
    Helper::halt( "cannot specify both keep and req" );

  bool req = param.has( "req" ) ;
  if ( param.has( "req" ) ) keeps = param.strset( "req" );
  
  if ( param.has( "drop" ) ) drops = param.strset( "drop" );
  
  if ( param.has( "keep" ) && param.has( "drop" ) )
    Helper::halt( "can only specify keep or drop with SIGNALS" );
  
  if ( ! ( param.has( "keep" ) || param.has( "drop" ) ) ) 
    Helper::halt( "need to specify keep or drop with SIGNALS" );

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
		Helper::halt( "could not find requested keep signal: " + *ss );
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


// CANONICAL

void proc_canonical( edf_t & edf , param_t & param )
{

  // just try to guess... no file specification
  if ( param.has( "guess" ) )
    {
      edf.guess_canonicals( param );
      return;
    }
  
  std::string file = param.requires( "file" );
  std::string group = param.requires( "group" );

  if ( ! param.has( "cs" ) )    
    edf.make_canonicals( file, group );
  else
    {
      const std::set<std::string> cs = param.strset( "cs" );
      edf.make_canonicals( file, group , &cs );
    }
}


// Adjust signals by other signals (similar to REFERENCE but different syntax)
void proc_adjust( edf_t & edf , param_t & param )
{
  edf.adjust( param );
}

// Reference tracks

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
  std::string new_channel = "";
  if ( make_new ) new_channel = param.value( "new" );
  int new_sr = 0;
  if ( make_new && param.has( "sr" ) ) new_sr = param.requires_int( "sr" );
  edf.reference( signals , references , make_new , new_channel , new_sr , false );
  
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
  std::string new_channel = "";
  if ( make_new ) new_channel = param.value( "new" );
  int new_sr = 0;
  if ( make_new && param.has( "sr" ) ) new_sr = param.requires_int( "sr" );
  
  edf.reference( signals , references , make_new , new_channel , new_sr , true );
  
}


// change record size for one or more signals

void proc_rerecord( edf_t & edf , param_t & param )
{
  double rs = param.requires_dbl( "dur" ); 

    logger << " altering record size from " << edf.header.record_duration << " to " <<  rs << " seconds\n";
  
    edf.reset_record_size( rs );

    logger << " now WRITE'ing EDF to disk, and will set 'problem' flag to skip to next EDF\n";

  proc_write( edf , param );
  globals::problem = true;
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
  edf.minmax( signals );
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





//
// Helper functions
//
	      
// void attach_annot( edf_t & edf , const std::string & astr )
// {
  
//   // are we checking whether to add this file or no? 
  
//   if ( globals::specified_annots.size() > 0 && 
//        globals::specified_annots.find( astr ) == globals::specified_annots.end() ) return;

//   // otherwise, annotation is either 
  
//   // 1) in memory (and may or may not be in a file,
//   // i.e. just created)
  
//   annot_t * a = edf.timeline.annotations( astr );
  
//   // 2) or, not in memory but can be retrieved from a file
  
//   if ( a == NULL )
//     {
      
//       // search for file for any 'required' annotation
//       // (e.g. annot=stages,spindles)
      
//       std::string annot_file = edf.annotation_file( astr );
      
//       if ( annot_file == "" ) 
// 	{
// 	  logger << " no instances of annotation [" 
// 		 << astr << "] for " << edf.id << "\n";		  
// 	}
//       else
// 	{
// 	  // add to list and load in data (note, if XML, load() does nothing
// 	  // as all XML annotations are always added earlier
	  	  
// 	  bool okay = annot_t::load( annot_file , edf );

// 	  if ( ! okay ) Helper::halt( "problem loading " + annot_file );
	  		  
// 	}
//     }
// }


void cmd_t::parse_special( const std::string & tok0 , const std::string & tok1 )
{


  // no console logging?
  if ( Helper::iequals( tok0 , "silent" ) ) 
    {
      globals::silent = Helper::yesno( tok1 );
      return;      
    }

  if ( Helper::iequals( tok0 , "verbose" ) )
    {
      globals::verbose = Helper::yesno( tok1 );
      return;
    }

  // specify indiv (i.e. can be used if ID is numeric)
  if ( Helper::iequals( tok0 , "id" ) )
    {
      globals::sample_list_id = tok1;
      return;
    }
  
  // dp for time output
  if ( Helper::iequals( tok0, "sec-dp" ) )
    {
      if ( ! Helper::str2int( tok1 , &globals::time_format_dp ) )
	Helper::halt( "expecting integer for sec-dp=N" );
      return;
    }
  
  // add signals?
  if ( Helper::iequals( tok0 , "sig" ) )
    {		  
      std::vector<std::string> tok2 = Helper::quoted_parse( tok1 , "," );		        
      for (int s=0;s<tok2.size();s++) 
	cmd_t::signallist.insert(Helper::unquote(tok2[s]));		  
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
    }

  // keep spaces
  if ( Helper::iequals( tok0 , "keep-spaces" ) )
    {
      globals::replace_channel_spaces = false;
      globals::replace_annot_spaces = false;
    }

  // keep spaces (annots only) 
  if ( Helper::iequals( tok0 , "keep-annot-spaces" ) )
    {
      globals::replace_annot_spaces = false;
    }
  
  // keep spaces (annots only) 
  if ( Helper::iequals( tok0 , "keep-channel-spaces" ) )
    {
      globals::replace_channel_spaces = false;
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
      globals::cmddefs.all_compressed( yesno );
      globals::cmddefs.none_compressed( !yesno );
      return;
    }

  // NSRR remapping
  if ( Helper::iequals( tok0 , "nsrr-remap" ) )
    {
      // clear pre-populated NSRR remapping
      if ( ! Helper::yesno( tok1 ) )
	nsrr_t::clear();
      return;
    }
  
  // generic annotation re-labelling, same format as 'alias'
  else if ( Helper::iequals( tok0 , "remap" ) )
    {
      nsrr_t::annot_remapping( tok1 );
      return;
    }

  // fix delimiter to tab only for .annot
  // default T --> tab-only=F is option to allow spaces  
  else if ( Helper::iequals( tok0 , "tab-only" ) )
    {
      globals::allow_space_delim = ! Helper::yesno( tok1 );
    }

  // default annot folder
  else if ( Helper::iequals( tok0 , "annot-folder" ) ||
	    Helper::iequals( tok0 , "annots-folder" ) ) 
    {
      if ( tok1[ tok1.size() - 1 ] != globals::folder_delimiter )
	globals::annot_folder = tok1 + globals::folder_delimiter ;
      else
	globals::annot_folder = tok1;		      
      return;
    }

  // if annot INST ID black, add hh:mm:ss
  else if ( Helper::iequals( tok0 , "inst-hms" ) )
    {
      globals::set_annot_inst2hms = Helper::yesno( tok1 );
      return;
    }

  // set INST ID to hh:mm:ss, whether it is blank or not
  else if ( Helper::iequals( tok0 , "force-inst-hms" ) )
    {
      globals::set_annot_inst2hms_force = Helper::yesno( tok1 );
      return;
    }

  // not enforce epoch check for .eannot
  else if ( Helper::iequals( tok0 , "no-epoch-check" ) )
    {
      globals::enforce_epoch_check = false; 
      return;
    }

  // set default epoch length
  else if ( Helper::iequals( tok0 , "epoch-len" ) )
    {
      if ( ! Helper::str2int( tok1 , &globals::default_epoch_len ) )
	Helper::halt( "epoch-len requires integer value, e.g. epoch-len=10" );
      return;
    }


  // additional annot files to add from the command line
  // i.e. so we don't have to edit the sample-list
  else if ( Helper::iequals( tok0 , "annots-file" ) ||
	    Helper::iequals( tok0 , "annots-files" ) ||
	    Helper::iequals( tok0 , "annot-file" ) ||
	    Helper::iequals( tok0 , "annot-files" ) )
    {
      globals::annot_files = Helper::parse( tok1 , "," );
      return;
    }


  // specified annots (only load these)
  else if ( Helper::iequals( tok0 , "annots" ) || Helper::iequals( tok0 , "annot" ) ) 
    {
      param_t dummy;     
      dummy.add( "dummy" , tok1 );
      globals::specified_annots = dummy.strset( "dummy" , "," );      
      return;
    }
  
  // signal alias?
  if ( Helper::iequals( tok0 , "alias" ) )
    {
      cmd_t::signal_alias( tok1 );
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

  // skip anyt EDF Annotations from EDF+
  if ( Helper::iequals( tok0 , "skip-edf-annots" ) )
    {
      globals::skip_edf_annots = Helper::yesno( tok1 );
      return;
    }

  // do not read ANNOT annotations
  if ( Helper::iequals( tok0 , "skip-annots" ) )
    {
      globals::skip_nonedf_annots = Helper::yesno( tok1 );
      return;
    }

  // do not read EDF or ANNOT annotations
  if ( Helper::iequals( tok0 , "skip-all-annots" ) )
    {
      globals::skip_edf_annots = globals::skip_nonedf_annots = Helper::yesno( tok1 );
      return;
    }

  // do not read FTR files 
  if ( Helper::iequals( tok0 , "ftr" ) )
    {
      globals::read_ftr = Helper::yesno( tok1 );
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

  
  // else a standard variable

  cmd_t::vars[ tok0 ] = tok1;

  
}



void cmd_t::define_channel_type_variables( edf_t & edf )
{
  // add these even if blank, so that scripts can always use
  
  std::string eeg = globals::list_channels( EEG , edf.header.label );
  //if ( eeg != "" )
  cmd_t::ivars[ edf.id ][ "eeg" ] = eeg;
    
  std::string ref = globals::list_channels( REF , edf.header.label );
  //if ( eeg != "" )
  cmd_t::ivars[ edf.id ][ "ref" ] = ref;

  std::string eog = globals::list_channels( EOG , edf.header.label );
  //if ( eog != "" )
  cmd_t::ivars[ edf.id ][ "eog" ] = eog;

  std::string ecg = globals::list_channels( ECG , edf.header.label );
  //if ( ecg != "" )
  cmd_t::ivars[ edf.id ][ "ecg" ] = ecg;

  std::string emg = globals::list_channels( EMG , edf.header.label );
  //if ( emg != "" )
    cmd_t::ivars[ edf.id ][ "emg" ] = emg;

  std::string leg = globals::list_channels( LEG , edf.header.label );
  //if ( leg != "" )
    cmd_t::ivars[ edf.id ][ "leg" ] = leg;

  std::string generic = globals::list_channels( GENERIC , edf.header.label );
  //if ( generic != "" )
    cmd_t::ivars[ edf.id ][ "generic" ] = generic;
  
  std::string airflow = globals::list_channels( AIRFLOW , edf.header.label );
  //if ( airflow != "" )
    cmd_t::ivars[ edf.id ][ "airflow" ] = airflow;

  std::string effort = globals::list_channels( EFFORT , edf.header.label );
  //if ( effort != "" )
    cmd_t::ivars[ edf.id ][ "effort" ] = effort;

  std::string oxygen = globals::list_channels( OXYGEN , edf.header.label );
  //if ( oxygen != "" )
    cmd_t::ivars[ edf.id ][ "oxygen" ] = oxygen;

  std::string position = globals::list_channels( POSITION , edf.header.label );
  //if ( position != "" )
    cmd_t::ivars[ edf.id ][ "position" ] = position;

  std::string light = globals::list_channels( LIGHT , edf.header.label );
  //if ( light != "" )
    cmd_t::ivars[ edf.id ][ "light" ] = light;

  std::string snore = globals::list_channels( SNORE , edf.header.label );
  //if ( snore != "" )
    cmd_t::ivars[ edf.id ][ "snore" ] = snore;

  std::string hr = globals::list_channels( HR , edf.header.label );
  // if ( hr != "" )
    cmd_t::ivars[ edf.id ][ "hr" ] = hr;

  std::string ignore = globals::list_channels( IGNORE , edf.header.label );
  //if ( ignore != "" )
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
// Attach i-vars from a file
//

void cmd_t::attach_ivars( const std::string & file )
{

  std::vector<std::string> files = Helper::parse( file , "," );

  for ( int f = 0; f < files.size() ; f++ )
    {
      std::string filename = Helper::expand( files[f] );

      if ( ! Helper::fileExists( filename ) )
	Helper::halt( "could not find " + filename );
      
      //std::cout << "reading from " << filename << "\n";
      
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
		Helper::halt( "inconsistent number of columns in " + filename );
	      
	      for (int c=0;c<ncols;c++)
		{
		  if ( c == idcol ) continue;		  
		  cmd_t::ivars[ tok[ idcol ] ][ head[ c ] ] = tok[c] ;		  
		  //std::cout << "setting: " << tok[ idcol ] << " " << head[ c ]  << "  " << tok[c]  << "\n";
		  
		}
	      
	    }

	  // read header now
	  header = false; 
	}
      IN1.close();
          
      //      logger << "  attached " << ncols - 1 << " from " << filename << "\n";
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
  specials.insert( "sec-dp" );
  specials.insert( "sig" ) ;
  specials.insert( "vars" );
  specials.insert( "ids" );    
  specials.insert( "add" ) ;
  specials.insert( "fail-list" ) ;
  specials.insert( "compressed" ) ;
  specials.insert( "nsrr-remap" ) ;
  specials.insert( "remap" ) ;
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
  specials.insert( "ftr" ) ;
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
