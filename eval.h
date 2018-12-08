
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

#ifndef __LUNA_EVAL_H__
#define __LUNA_EVAL_H__

#include <string>
#include <map>
#include <set>
#include <vector>
#include <sstream>
#include "helper/helper.h"

struct edf_t;

void attach_annot( edf_t & edf , const std::string & astr );

//
// Helper to parse command syntax
//

class param_t
{

 public:
  
  void add( const std::string & option , const std::string & value = "" ) { opt[ option ] = value; }  
  
  void parse( const std::string & s )
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
  
  void update( const std::string & id , const std::string & wc )
  {
    // replace all instances of 'globals::indiv_wildcard' with 'id'
    // for all values
    std::map<std::string,std::string>::iterator ii = opt.begin();
    while ( ii != opt.end() ) 
      {
	std::string v = ii->second;
	bool changed = false;
	while ( v.find( wc ) != std::string::npos )
	  {
	    int p = v.find( wc );
	    v = v.substr( 0 , p ) + id + v.substr(p+1);
	    changed = true;
	  }
	if ( changed ) ii->second = v;
	++ii;
      }

  }

  void clear() { opt.clear(); } 
  
  bool has(const std::string & s ) const { return opt.find(s) != opt.end(); } 

  std::string match( const std::string & s ) 
    {
      std::string ms = "";
      int m = 0;
      std::map<std::string,std::string>::iterator ii = opt.begin();
      while ( ii != opt.end() )
	{
	  if ( ii->first.substr( 0, s.size() ) == s ) 
	    {
	      ms = ii->first;
	      ++m;
	    }
	  ++ii;
	}
      if ( m == 1 ) return ms;
      return "";
    }
  
  std::string value( const std::string & s ) const { return has(s) ? opt.find(s)->second : "" ; }
  
  bool single() const { return opt.size() == 1; }
  
  std::string single_value() const 
    { if ( ! single() ) Helper::halt( "no single value" ); return opt.begin()->first; }

  std::string requires( const std::string & s ) const
    {
      if ( ! has(s) ) Helper::halt( "command requires parameter " + s );
      return value(s);
    }
  
  int requires_int( const std::string & s ) const
  {
    if ( ! has(s) ) Helper::halt( "command requires parameter " + s );
    int r;
    if ( ! Helper::str2int( value(s) , &r ) ) 
      Helper::halt( "command requires parameter " + s + " to have an integer value" );
    return r;
  }
  
  double requires_dbl( const std::string & s ) const
  {
    if ( ! has(s) ) Helper::halt( "command requires parameter " + s );
    double r;
    if ( ! Helper::str2dbl( value(s) , &r ) ) 
      Helper::halt( "command requires parameter " + s + " to have a numeric value" );
    return r;
  }

  std::string dump( const std::string & indent = "  ", const std::string & delim = "\n" ) const
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

  std::set<std::string> strset( const std::string & k , const std::string delim = "," ) const
    {
      std::set<std::string> s;
      if ( ! has(k) ) return s;
      std::vector<std::string> tok = Helper::quoted_parse( value(k) , delim );
      for (int i=0;i<tok.size();i++) s.insert( Helper::unquote( tok[i]) );
      return s;
    }
  
  std::vector<std::string> strvector( const std::string & k , const std::string delim = "," ) const
    {
      std::vector<std::string> s;
      if ( ! has(k) ) return s;
      std::vector<std::string> tok = Helper::quoted_parse( value(k) , delim );
      for (int i=0;i<tok.size();i++) s.push_back( Helper::unquote( tok[i]) );
      return s;
    }

  std::vector<double> dblvector( const std::string & k , const std::string delim = "," ) const
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

  std::vector<int> intvector( const std::string & k , const std::string delim = "," ) const
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

  std::set<std::string> keys() const
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

 private:

  std::map<std::string,std::string> opt;

};

class cmd_t 
{
  
 public:
  
  cmd_t() 
    {
      reset();
      error = ! read();
    }

  cmd_t( const std::string & str ) 
    {
      reset();
      error = ! read( &str , true ); 
    }

  static void add_cmdline_cmd( const std::string & c ) 
  {
    cmdline_cmds.append( c + " " );
  }
  
  // by default (str == NULL) read from command line arguments, 
  // then from STDIN;  otherwise, parse the commands in the 
  // text string 

  bool read( const std::string * str = NULL , 
	     bool silent = false ); 
  
  void reset() 
  {
    cmds.clear();
    params.clear();
    line = "";
    error = false;
    will_quit = false;
  }

  void replace_wildcards( const std::string & id );

  bool eval( edf_t & ) ;
  
  bool empty() const { return will_quit; }

  bool valid() const 
  {    
    if ( error ) return false;
    /* for (int c=0;c<cmds.size();c++) */
    /*   if ( commands.find( cmds[c] ) == commands.end() ) return false; */
    return true;
  }

  bool badline() const { return error; } 

  std::string offending() const { return ( error ? line : "" ); }
  
  int num_cmds() const { return cmds.size(); }
  
  std::string cmd(const int i) { return cmds[i]; }
    
  param_t & param(const int i) { return params[i]; }
  
  bool process_edfs() const
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
  
  bool is( const int n , const std::string & s ) const
  {
    if ( n < 0 || n >= cmds.size() ) Helper::halt( "bad command number" );
    return Helper::iequals( cmds[n] , s );
  }
  
  std::string data() const { return input; } 

  bool quit() const { return will_quit; }
  void quit(bool b) { will_quit = b; }

  //
  // Static members (i.e. from command line)
  //

  static std::string                        input;  
  static std::map<std::string,std::string>  vars;

  static std::string                        cmdline_cmds;
  static std::string                        stout_file;
  static bool                               append_stout_file;

  // command-specific parameters (i.e. from command-file)
  static std::set<std::string>    signallist;
  
  static std::map<std::string,std::string> label_aliases;
  static std::map<std::string,std::vector<std::string> > primary_alias;

  static void signal_alias( const std::string & s )
  {
    // format canonical|alias1|alias2 , etc.
    std::vector<std::string> tok = Helper::quoted_parse( s , "|" );    
    if ( tok.size() < 2 ) Helper::halt( "bad format for signal alias:  canonical|alias 1|alias 2" );
    const std::string primary = Helper::unquote( tok[0] );
    for (int j=1;j<tok.size();j++) 
      {
	label_aliases[ Helper::unquote( tok[j] ) ] = primary;
	primary_alias[ primary ].push_back( Helper::unquote( tok[j] ) );
      }
    
  }


  static const std::set<std::string> & signals() { return signallist; }
  static void clear_signals() { signallist.clear(); }

  static std::string signal_string() 
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


 private:
  
  // a command is a single input (0+ EDFs) linked to 1 or more
  // commands, that are performed sequentially on each EDF

  std::string line; // raw input from STDIN, stored to report in case of error

  bool error;
  bool will_quit;
  
    
  std::vector<std::string> cmds;
  std::vector<param_t>     params;
  std::vector<param_t>     original_params;
  
  static std::set<std::string> commands;
  
  static void populate_commands() 
  {
    commands.insert( "VALIDATE" );
    commands.insert( "HEADERS" );
    commands.insert( "SUMMARY" );
    commands.insert( "PSD" );
    commands.insert( "SPINDLES" );
    commands.insert( "RMS" );
    commands.insert( "MASK" );
    commands.insert( "CWT" );
    commands.insert( "RESCALE" );
    commands.insert( "DUMP" );
    commands.insert( "DUMP-MASK" );
    commands.insert( "DUMMY" );
    commands.insert( "ARTIFACTS" );
    commands.insert( "TIME-TRACK" );
  }

};



//
// Primary commands
//

void proc_summaries( const std::string & , const std::string & , int , cmd_t & , param_t & param );
void proc_summaries( edf_t & , param_t & );
void proc_desc( edf_t & , param_t & );

void proc_validate( edf_t & , param_t & );
void proc_list_annots( const std::string & edffile , const std::string & rootname , const std::vector<std::string> & );
void proc_list_all_annots( edf_t & , param_t & );



void proc_tag( param_t & );
void set_tag( const std::string & t = "." );
void proc_anon( edf_t & , param_t & );
void proc_write( edf_t & , param_t & );
void proc_restructure( edf_t & , param_t & );
void proc_drop_signals( edf_t & , param_t & );
void proc_scale( edf_t & , param_t & , const std::string & s );
void proc_flip( edf_t & , param_t & );
void proc_reference( edf_t & , param_t & );
void proc_rerecord( edf_t & edf , param_t & param );

void proc_dump( edf_t & , param_t & );
void proc_dump_mask( edf_t & , param_t & );
void proc_file_mask( edf_t & , param_t & );
void proc_file_annot( edf_t & , param_t & );
void proc_sleep_stage( edf_t & , param_t & );

void proc_record_dump( edf_t & , param_t & );
void proc_intervals( param_t & , const std::string & );
void proc_epoch_dump( edf_t & , param_t & );
void proc_epoch_matrix( edf_t & , param_t & );
void proc_epoch_mask( edf_t & , param_t & );
void proc_mask( edf_t & , param_t & );
void proc_epoch( edf_t & , param_t & );
void proc_slice( edf_t & , param_t & , int );

void proc_timetrack( edf_t & , param_t & );

void proc_covar( edf_t & , param_t & );

void proc_artifacts( edf_t & , param_t & );
void proc_rms( edf_t & , param_t & );
void proc_mse( edf_t & , param_t & );
void proc_lzw( edf_t & , param_t & );
void proc_zratio( edf_t & , param_t & );

void proc_resample( edf_t & , param_t & );
void proc_filter( edf_t & , param_t & );
void proc_filter_legacy( edf_t & , param_t & );
void proc_filter_design( edf_t & , param_t & );
void proc_filter_design_cmdline();
void proc_cwt_design_cmdline();
void proc_psd( edf_t & , param_t & );
void proc_1overf_norm( edf_t & , param_t & );
void proc_tv_denoise( edf_t & , param_t & );

void proc_fiplot( edf_t & , param_t & );
void proc_mi( edf_t & , param_t & );
void proc_ica( edf_t & , param_t & );
void proc_emd( edf_t & , param_t & );
void proc_leave_one_out( edf_t & , param_t & );
void proc_correl( edf_t & , param_t & );
void proc_coh( edf_t & , param_t & );
void proc_coh_legacy( edf_t & , param_t & );
void proc_elec_distance( edf_t & , param_t & );

void proc_pac( edf_t & , param_t & );
void proc_cfc( edf_t & , param_t & );

void proc_ecgsuppression( edf_t & , param_t & );
void proc_bpm( edf_t & , param_t & );

void proc_spindles( edf_t & , param_t & );
void proc_polarity( edf_t & , param_t & );
void proc_slowwaves( edf_t & , param_t & );
void proc_cwt( edf_t & , param_t & );
void proc_spike( edf_t & , param_t & );



#endif
