

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
#include <vector>
#include <set>

struct param_t;


//
// cmddefs_t handles a) help, listing parametes and output for each command, 
//                   b) groups commands by domains
//                   c) provides a validation of known parameters for a command check()
//                   d) tracks whether a given output table (or all) should be dumped as a plain-text file or no  out_plaintext()


struct strata_t;

struct tfac_t { 

  // need to handle: baseline strata (no factors),    augemnted TAGs (i.e something plus TAGs) 
  // still has same rules
  // RULE::: TAGs will always start

  // The following are *ignored* when constructing a tfac_t: 
  //   factors that start with underscore '_', as these as command names
  //   factors that are listed as TAGs here
  // i.e. as we already know command name, and the TAGs will not define the variables to use
  
  tfac_t( const std::string & s , const std::string & delim = "," ) ;
  
  tfac_t( const strata_t & s );

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

#ifdef __clang__
  void __attribute__((optnone)) init();
#else
#ifdef __GNUC__    
  void __attribute__((optimize(0))) init();
#else
  void init();
#endif
#endif


  
  // domain description
  void add_domain( const std::string & domain , const std::string & label ,  const std::string & desc );

  bool is_domain( const std::string & d ) ;

  // command description 
  void add_cmd( const std::string & domain , const std::string & cmd , const std::string & desc , const bool hide = false );


  // hidden command description 
  void hidden_cmd( const std::string & domain , const std::string & cmd , const std::string & desc );
  
  bool is_cmd( const std::string & c ) ;
  
  // command URLs , e.g.  zzz.bwh.harvard.edu/luna/ref/
  void add_url( const std::string & cmd , const std::string & url ) ;

  void add_note( const std::string & cmd , const std::string & note ) ;

  // parameters for this command
  void add_param( const std::string & cmd , const std::string & param , 
		  const std::string & ex ,  // "" if none
  		  const std::string & desc , 
		  const std::string & requirements = "" ,
		  const bool hide = false );

  // hide parameter for this command
  void hidden_param( const std::string & cmd , const std::string & param , 
		   const std::string & ex ,  // "" if none
		   const std::string & desc , 
		   const std::string & requirements = "" );
  
  // output from this command , "CMD" , "F,B,CH,E" , "desc" , is compressed Y/N
  void add_table( const std::string & cmd , const std::string & factors , const std::string & desc , bool isz = false , bool hide = false );

  void hidden_table( const std::string & cmd , const std::string & factors , const std::string & desc , bool isz = false );

  // ensure table (from REPORT)
  void ensure_table( const std::string & cmd , const std::string & factors );

  // add variable
  void add_var( const std::string & cmd , const std::string & factors , const std::string & var , const std::string & desc , const bool hide = false );

  // add hidden variable
  void hidden_var( const std::string & cmd , const std::string & factors , const std::string & var , const std::string & desc );

  // register extra col for output (controlled by caller)
  void register_var( const std::string & cmd , const std::string & factors , const std::string & var , const bool value = true );
  
  
  //
  // check parameters
  //

  bool check( const std::string & cmd ) const ; // command
  
  bool check( const std::string & cmd , const std::set<std::string> & k , std::set<std::string> * ) const ; // params


  //
  // add 'hidden' commands cmd/tables/vars/param
  //

  void add_cmd1( const std::string & domain ,
		 const std::string & cmd ,
		 const std::string & desc ,
		 const bool hide = false );

  
  void add_param1( const std::string & cmd , const std::string & param , 
		   const std::string & ex ,  // "" if none
		   const std::string & desc , 
		   const std::string & requirements = "" ,
		   const bool hide = false );

  
  void add_table1( const std::string & cmd ,
		   const std::string & factors ,
		   const std::string & desc ,
		   bool isz = false ,
		   bool hide = false );
  
  void add_var1( const std::string & cmd ,
		 const std::string & factors ,
		 const std::string & var ,
		 const std::string & desc ,
		 const bool hide = false );

  
  //
  // show help
  //

  // describe one command 
  std::string help( const std::string & cmd , bool show_domain_label = true , bool verbose = false , bool primary = false ) const;

  // list all domains
  std::string help_domains() const;

  // desc for a domain
  std::string help_domain( const std::string & d ) const;
  
  // list all commands in a domain
  std::string help_commands( const std::string & d , const bool primary ) const;

  // list all commands 
  std::string help_commands() const;

  //
  // programmatically report all commands/options
  // 

  // lists: all domains
  //      : all commands in a given domain
  //      : all param for a given command
  //      : all tables for a given command
  //      : all vars for a given command/table

  // queries (desc)
  //      : for a given domain
  //      : for a given command
  //      : for a given command/param
  //      : for a given command/table
  //      : for a given command/table/var
  
  std::vector<std::string> fetch_doms( const bool all = true ) const;
  std::vector<std::string> fetch_cmds( const std::string & dom , const bool all = true ) const;
  std::vector<std::string> fetch_params( const std::string & cmd, const bool all = true ) const;
  std::vector<std::string> fetch_tbls( const std::string & cmd, const bool all = true ) const;
  std::vector<std::string> fetch_vars( const std::string & cmd, const std::string & tbl, const bool all = true ) const;

  std::string fetch_desc_dom( const std::string & dom ) const;
  std::string fetch_desc_cmd( const std::string & cmd ) const;
  std::string fetch_desc_param( const std::string & cmd, const std::string & param ) const;
  std::string fetch_desc_tbl( const std::string & cmd, const std::string & tbl ) const;
  std::string fetch_desc_var( const std::string & cmd, const std::string & tbl, const std::string & var ) const;

  
  //
  // hide/show variables
  //

  void hide_all();
  void hide_cmd( const std::string & cmd );

  void hide_table( const std::string & cmd , const std::string & factors );
  void hide_table( const std::string & cmd , const tfac_t & factors );

  void hide_var( const std::string & cmd , const std::string & factors , const std::string & var );
  void hide_var( const std::string & cmd , const tfac_t & factors , const std::string & var );

  void show_all( const bool status = true );
  void show_cmd( const std::string & cmd , const bool status = true );

  void show_table( const std::string & cmd , const std::string & factors , const bool status = true );
  void show_table( const std::string & cmd , const tfac_t & factors , const bool status = true );

  void show_var( const std::string & cmd, const std::string & factors, const std::string & var, const bool status = true );
  void show_var( const std::string & cmd, const tfac_t & factors, const std::string & var, const bool status = true );
  
  

  
  //
  // indicate whether table should be compressed or no
  //

  void all_compressed( bool b ) ;
  bool all_compressed() const ;

  void none_compressed( bool b );
  bool none_compressed() const ;


  // allow change after table has been registered  
  void set_compressed( const std::string & cmd , 
		       const tfac_t & tfac , 
		       const bool b = true
		       );

  void set_compressed( const std::string & cmd , 		       
		       const bool b = true
		       );
  
  bool out_compressed( const std::string & cmd , 
		       const tfac_t & tfac ) const;
  
  

  bool exists( const std::string & cmd ,
	       const tfac_t & tfac ) const;

  
  //
  // list factors
  //
  
  void add_tag( const std::string & tag );
  void clear_tags();
  bool is_tag( const std::string & tag ) const ;
  std::set<std::string> variables( const std::string & cmd , const param_t * param , const tfac_t & tfac  );

 private:

  
  // domain->human label    ->desc
  std::map<std::string,std::string> domain_label;
  std::map<std::string,std::string> domain_desc;
  std::vector<std::string> domain_ordered;
  
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
  // primary cmd/tables/params/vars (i.e. for help function)
  //
  
  bool is_primary_cmd( const std::string & cmd ) const;
  bool is_primary_par( const std::string & cmd , const std::string & param ) const;
  bool is_primary_tbl( const std::string & cmd , const tfac_t & tfac ) const;
  bool is_primary_var( const std::string & cmd , const tfac_t & tfac , const std::string & var ) const;

  std::set<std::string> pri_cmd;
  std::map<std::string,std::set<std::string> > pri_par;
  std::map<std::string,std::set<tfac_t> > pri_tbl;
  std::map<std::string,std::map<tfac_t,std::set<std::string> > > pri_var;
  
  //
  // output
  //
	      
  // cmd->table->desc
  std::map<std::string,std::map<tfac_t,std::string> > otables;
  
  // cmd->{table-set}->compressed table T/F
  std::map<std::string,std::map<tfac_t,bool> > ofacs;
  
  // cmd->table->var->desc
  std::map<std::string,std::map<tfac_t,std::map<std::string,std::string> > > ovars;

  // cmd->table->var->tout
  std::map<std::string,std::map<tfac_t,std::map<std::string,std::string> > > otout;

  //
  // hidden status (i.e. not reported in output)
  //

  std::map<std::string,bool> chide;  // cmds
  std::map<std::string,std::map<std::string,bool> > phide;  // parameters
  std::map<std::string,std::map<tfac_t,bool> > ohide;  // tables
  std::map<std::string,std::map<tfac_t,std::map<std::string,bool> > > vhide;  // variables
      
  bool is_hidden_cmd( const std::string & c ) const;
  
  bool is_hidden_param( const std::string & c , const std::string & p ) const;

  bool is_hidden_table( const std::string & c , const tfac_t & tfac ) const;

public:
  bool is_hidden_var( const std::string & c , const tfac_t & tfac , const std::string & v ) const;
private:
  
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


