

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

#ifndef __LUNA_CMDDEFS_H__
#define __LUNA_CMDDEFS_H__

#include <string>
#include <map>
#include "eval.h"
#include "helper/helper.h"
#include <vector>
#include <set>
#include <sstream>



//
// cmddefs_t handles a) help, listing parametes and output for each command, 
//                   b) groups commands by domains
//                   c) provides a validation of known parameters for a command check()
//                   d) tracks whether a given output table (or all) should be dumped as a plain-text file or no  out_plaintext()


struct tfac_t { 

  // need to handle: baseline strata (no factors),    augemnted TAGs (i.e something plus TAGs) 
  // still has same rules
  // RULE::: TAGs will always start

  // The following are *ignored* when constructing a tfac_t: 
  //   factors that start with underscore '_', as these as command names
  //   factors that are listed as TAGs here
  // i.e. as we already know command name, and the TAGs will not define the variables to use
  
  tfac_t( const std::string & s ) ;

  std::string as_string( const std::string & delim = "," ) const;

  // LT operator
  bool operator< ( const tfac_t & rhs ) const;

  // EQ operator
  bool operator== ( const tfac_t & rhs ) const;

  // key (non-command, non-tag) factors that represent the virtual table
  std::set<std::string> fac;

};



class cmddefs_t 
{
  
 public:
  
  cmddefs_t();
  
  void init();

  // domain description
  void add_domain( const std::string & domain , const std::string & label ,  const std::string & desc )
  {
    domain_label[ domain ] = label;
    domain_desc[ domain ] = desc;
  }

  bool is_domain( const std::string & d ) 
  {
    return domain_label.find( d ) != domain_label.end();
  }

  // command description 
  void add_cmd( const std::string & domain , const std::string & cmd , const std::string & desc )
  {
    dcmds[ domain ].insert( cmd );
    cmds[ cmd ] = desc ; 
    cdomain[ cmd ] = domain;
  }

  bool is_cmd( const std::string & c ) 
  {
    return cmds.find( c ) != cmds.end();
  }
  
  
  // command URLs , e.g.  zzz.bwh.harvard.edu/luna/ref/
  void add_url( const std::string & cmd , const std::string & url ) 
  {
    if ( cmds.find( cmd ) == cmds.end() ) Helper::halt( cmd + " not registered" );
    curl[ cmd ] = url;
  }

  void add_note( const std::string & cmd , const std::string & note ) 
  {
    if ( cmds.find( cmd ) == cmds.end() ) Helper::halt( cmd + " not registered" );
    cnotes[ cmd ] = note;
  }

  // parameters for this command
  void add_param( const std::string & cmd , const std::string & param , 
		  const std::string & ex ,  // "" if none
  		  const std::string & desc , 
		  const std::string & requirements = "" )
  {
    pdesc[ cmd ][ param ] = desc;
    preq[ cmd ][ param ] = requirements;
    px[ cmd ][ param ] = ex;
  }
  
  // output from this command , "CMD" , "F,B,CH,E" , "desc" , is compressed Y/N
  void add_table( const std::string & cmd , const std::string & factors , const std::string & desc , bool isz = false )
  {
    tfac_t tfac( factors );
    otables[ cmd ][ tfac ] = desc ; 
    ofacs[ cmd ][ tfac ] = isz ; 
  }
  
  // add variable
  void add_var( const std::string & cmd , const std::string & factors , const std::string & var , const std::string & desc )
  {
    tfac_t tfac( factors );
    ovars[ cmd ][ tfac ][ var ] = desc;
  }

  
  //
  // check parameters
  //

  bool check( const std::string & cmd ) const ; // command
  
  bool check( const std::string & cmd , const std::set<std::string> & k , std::set<std::string> * ) const ; // params

  //
  // show help
  //

  // describe one command 
  std::string help( const std::string & cmd , bool show_domain_label = true , bool verbose = false ) const;

  // list all domains
  std::string help_domains() const;
  
  // list all commands in a domain
  std::string help_commands( const std::string & d ) const;

  // list all commands 
  std::string help_commands() const;


  //
  // indicate whether table should be compressed or no
  //

  void all_compressed( bool b ) { allz = b; } 
  bool all_compressed() const { return allz; }

  void none_compressed( bool b ) { nonez = b; } 
  bool none_compressed() const { return nonez; }


  // allow change after table has been registered  
  void set_compressed( const std::string & cmd , 
		       const tfac_t & tfac , 
		       const bool b = true
		       );
  
  bool out_compressed( const std::string & cmd , 
		       const tfac_t & tfac ) const;
  
  

  bool exists( const std::string & cmd ,
	       const tfac_t & tfac ) const;

  
  //
  // list factors
  //
  
  void add_tag( const std::string & tag ) { tags.insert( tag ); } 
  void clear_tags() { tags.clear(); }
  bool is_tag( const std::string & tag ) const { return tags.find( tag ) != tags.end(); } 
  std::set<std::string> variables( const std::string & cmd , const param_t * param , const tfac_t & tfac  );

 private:

  
  // domain->human label    ->desc
  std::map<std::string,std::string> domain_label;
  std::map<std::string,std::string> domain_desc;
  
  // domain->cmds
  std::map<std::string,std::set<std::string> > dcmds;  

  // cmd->desc
  std::map<std::string,std::string> cmds;

  // cmd->notes
  std::map<std::string,std::string> cnotes;

  // cmd->url
  std::map<std::string,std::string> curl;

  std::string url_root;
  
  // cmd->domain
  std::map<std::string,std::string> cdomain;
  
  //cmd->param->param->example
  std::map<std::string,std::map<std::string,std::string> > pdesc;

  //cmd->param->param->example
  std::map<std::string,std::map<std::string,std::string> > px;

  // cmd->param->param->requirements (free-text)
  std::map<std::string,std::map<std::string,std::string> > preq;

  //
  // output
  //
	      
  // cmd->table->desc
  std::map<std::string,std::map<tfac_t,std::string> > otables;
  
  // cmd->{table-set}->compressed table T/F
  std::map<std::string,std::map<tfac_t,bool> > ofacs;
  
  // cmd->table->var->desc
  std::map<std::string,std::map<tfac_t,std::map<std::string,std::string> > > ovars;
  
  // all, or no, output should be compressed
  bool allz;
  bool nonez;

  //
  // TAGs, need to keep track of these so they can be ignored when determining 
  // the table to write to (i.e. these do not determine which variables will be 
  // present, they will just add extra strata/rows
  //

  std::set<std::string> tags;
  
};


#endif


