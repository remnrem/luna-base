
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

#ifndef __MERGE_H__
#define __MERGE_H__

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <set>

#include "merge-helpers.h"

struct options_t {
  options_t() {
    skip_folders.insert( "extra" );

    missing_data_outsymbol = "NA";

    missing_data_symbol.insert(  "NA" );
    missing_data_symbol.insert(  "?" );
    missing_data_symbol.insert(  "." );
    missing_data_symbol.insert(  "NaN" );
    missing_data_symbol.insert(  "nan" );
    missing_data_symbol.insert(  "-nan" );

    domain_includes.clear();
    file_excludes.clear();
    var_excludes.clear();
    hms_delim = ":.";
    date_delim = "/-.";
    show_fac = true;
    verbose = false;
    strict = false;
    assume_txt = true;
    max_var_len = 100;
    numeric_strata_encoding = false;
  }

  bool verbose;
  bool strict; // F skip bad stuff / T halt if bad stuff found
  bool assume_txt; // only look at .txt files (i.e. requires that extension)
  int max_var_len;
  bool numeric_strata_encoding; // VAR.1, VAR.2 instead of VAR.FAC_LVL_FAC_LVL
  std::set<std::string> skip_folders;
  std::string missing_data_outsymbol; // NA for output
  std::set<std::string> missing_data_symbol; // for input
  std::map<std::string,std::set<std::string> > domain_includes;
  std::set<std::string> file_excludes;
  std::set<std::string> var_excludes;

  // hh:mm:ss delimiter characters (allowed up to two)
  std::string hms_delim;
  // date delimiter chars
  std::string date_delim;

  // show factors in varnames
  bool show_fac;
  
  bool is_missing( const std::string & val ) const
  {
    std::set<std::string>::const_iterator ii = missing_data_symbol.begin();
    while ( ii != missing_data_symbol.end() )
      {
	if ( iequals( *ii , val ) ) return true;
	++ii;
      }
    return false;
  }

  bool read_file( const std::string & filetag ) const
  {
    return file_excludes.find( filetag ) == file_excludes.end() ;
  }
  
  bool read_domain( const std::string & domain , const std::string & group ) const {

    if ( domain_includes.size() == 0 ) return true;
    std::map<std::string,std::set<std::string> >::const_iterator ii = domain_includes.find( domain );

    // domain not in include list
    if ( ii == domain_includes.end() ) return false;

    // domain in include list, and no groups specified
    if ( ii->second.size() == 0 ) return true;

    // groups specified, is the domain::group pairing in the include
    return ii->second.find( group ) != ii->second.end() ;

  }

  int parse_opt( const std::string & o , std::string * k , std::string * v )
  {
    if ( o.size() <= 1 ) return 0;
    std::vector<std::string> tok = parse( o.substr(1) , "=" );
    if ( tok.size() == 0 || tok.size() > 2 ) return 0;
    *k = tok[0];
    if ( tok.size() == 2 ) *v = tok[1];
    return tok.size();
  }

  void include_domain( const std::string & t )
  {
    std::vector<std::string> tok = parse( t , "_" );
	  
    if ( tok.size() > 2 || tok.size() == 0 ) 
      halt( "invalid domain_group specification: " + t );
	  
    if ( tok.size() == 1 ) 
      {
	std::set<std::string> empty;
	domain_includes[ tok[0] ] = empty;
      }
    else
      domain_includes[ tok[0] ].insert( tok[1] );
  }
  
};


extern options_t options;

enum type_t { ROWID , FACTOR , TEXT , INT , FLOAT , YESNO , DATE , TIME };

bool type_check( const std::string & value , type_t );


struct alias_t {

  void add_canonical( const std::string & canonical )
  {
    canonicals.insert( canonical );
  }

  void add_alias( const std::string & alias , const std::string & canonical ) 
  {
    // nonsensical
    if ( alias == canonical ) halt( "alias and canonical equal: " + alias );

    // check canonical term not already an alias
    if ( canonicals.find( alias ) != canonicals.end() )
      halt("cannot specify " + alias + " as both an ALIAS and canonical term" );

    // check alias not already pointing to a different canonical term
    std::map<std::string,std::string>::const_iterator aa = a.find( alias );
    if ( aa != a.end() && aa->second != canonical )
      halt( "alias " + alias + " cannot point to multiple canonical values (" + canonical + " and " + aa->second );

    a[ alias ] = canonical;
    canonicals.insert( canonical );

  }
  
  std::string unalias( const std::string & n ) const
  {
    std::map<std::string,std::string>::const_iterator aa = a.find(n);
    if ( aa == a.end() ) return n;
    return aa->second;
  }

  // alias --> canonical
  std::map<std::string,std::string> a;

  // check: canonical cannot also be an alias
  std::set<std::string> canonicals;
  
};

struct domain_t;

struct var_t {

  var_t() { }

  // for lookups
  var_t( const std::string & name ) : name(name) { } 
	   
  var_t( const domain_t & domain ,
	 const std::string & name ,
	 const std::string & t , 
	 const std::string & label ,
	 const std::string & b = "" );

  std::string print_type() const {    
    if ( type == FACTOR ) return "Factor";
    else if ( type == INT ) return "Integer";
    else if ( type == FLOAT ) return "Numeric";
    else if ( type == YESNO ) return "YesNo";
    else if ( type == DATE ) return "Date";
    else if ( type == TIME ) return "Time";
    else return "Text";
  }
  
  bool operator<( const var_t & rhs ) const {
    return name < rhs.name ; 
  }
  
  std::string name; // e.g. DENS.F_11.5_CH_C3_SS_N2
  std::string base; // e.g. DENS
  std::string label;
  type_t type;
  std::string domain_name, domain_group;
  std::map<std::string,std::string> fac2lvl; // track (for data dictionary output)
};


struct domain_t {

  domain_t( const std::string & name , const std::string & group )
    : name(name) , group(group) 
  {
    missing.clear();
  }
  
  domain_t( const std::string & filename )
  {
    read( filename );
  }
  
  int read( const std::string & filename );  

  std::string name; // e.g. 'eeg'
  std::string group; // e.g. 'spindles'
  // combined to eeg_spindles : i.e. enfore underscore delimiter

  std::map<std::string,var_t> variables;

  // domain-specific aliases for variable names
  alias_t aliases;
  
  // domain-specific missing data symbol(s)?
  std::set<std::string> missing;
  
  bool has( const std::string & varname) const
  {
    const var_t * var = variable( varname );
    if ( var == NULL ) return false;
    return true;
  }

  bool has( const std::string & varname , type_t type ) const
  {
    const var_t * var = variable( varname );
    if ( var == NULL ) return false;
    return var->type == type;
  }
    
  const var_t * variable( const std::string & varname ) const
  {
    std::map<std::string,var_t>::const_iterator ii = variables.find( varname );
    if ( ii == variables.end() ) return NULL;
    return &(ii->second);
  }
  
  bool operator<( const domain_t & rhs ) const {
    if ( name < rhs.name ) return true;
    if ( name > rhs.name ) return false;
    return group < rhs.group ; 
  }

};





struct value_t {
  value_t() { } 
  value_t( const std::string & data ) : data(data) { }  
  std::string data;
  // can add type checks here, and missing codes, etc
};



struct indiv_t {

  indiv_t( const std::string & id ) : id(id) { } 

  std::string print() const {
    std::stringstream ss;
    ss << id << "\n";
    std::map<var_t,value_t>::const_iterator ii = values.begin();
    while ( ii != values.end() )
      {
	ss << ii->first.name << " --> " << ii->second.data << "\n";
	++ii;
      }
    return ss.str();
  }
  
  // merge two datarows
  void add( const indiv_t & indiv ) {
    if ( id != indiv.id ) halt( "trying to merge different IDs" );
    std::map<var_t,value_t>::const_iterator ii = indiv.values.begin();
    while ( ii != indiv.values.end() )
      {
	if ( values.find( ii->first ) != values.end() )
	  halt( "multiple obervations for " + id
		+ " for variable: " + ii->first.name );
	values[ ii->first ] = ii->second;
	++ii;
      }
  }

  void add( const var_t & var , const value_t & value )
  {
    if ( values.find( var ) != values.end() )
      halt( "multiple values for " + id + " " + var.name );
    values[ var ] = value;
  }
  
  std::string id;
  
  std::map<var_t,value_t> values; 

  bool operator<( const indiv_t & rhs ) const {
    return id < rhs.id;
  }
  
};



struct dataset_t {

  void add( const domain_t & domain )
  {
    std::cerr << " ++ adding domain " << domain.name
	      << "::" << domain.group
	      << " (" << domain.variables.size() << " variables)\n";
    domains.insert( domain );
  }
  
  void read( const std::string & filename );
  
  void write( std::ofstream & );
  
  std::set<indiv_t> indivs;

  std::set<std::string> files;

  std::set<domain_t> domains;

  std::set<var_t> xvars; // expanded vars, e.g. DENS.F_11_CH_C3_SS_N2, or DENS.1, DENS.2 etc
  
  std::map<std::string,int> obscount; // count non-missing obs
  
  std::map<std::string,std::string> faclabels; // track labels for factors (for dict output)

  std::map<std::string,int> strata2number; // for numeric strata encoding
  
  var_t xvar( const var_t & var , const std::vector<std::string> & fac , const std::vector<std::string> & lvl )
  {

    if ( fac.size() != lvl.size() )
      halt( "internal error in add_xvar" );

    //
    // no expansion needed?
    //
    
    if ( fac.size() == 0 ) {
      xvars.insert( var );      
      return var;
    }

    //
    // already done expansion? HMM.. this is problematic, e.g. if the same ROOT variable
    // can feature with different strata (e.g. NREM_MINS baseline and NREM_MINS/C
    // so, skip this step, always force evaluation
    //

    if ( false ) 
      {
	if ( xvars.find( var ) != xvars.end() )
	  return var;
      }

    //
    // expand variable name given factors
    //

    std::string xname = var.name;
    std::string xlabel = var.label + " (";

    var_t xv;
    
    for (int i=0;i<fac.size();i++)
      {
	std::string fac2 = search_replace( fac[i] , '_' , '.' );
	std::string lvl2 = search_replace( lvl[i] , '_' , '.' );

	// FAC_LVL or just LVL encoding in varnames
	if ( options.show_fac )
	  xname += (i > 0 ? "_" : "." ) + fac2 + "_" + lvl2;
	else
	  xname += (i > 0 ? "_" : "." ) + lvl2;
	
	xlabel += (i > 0 ? ", " : "" ) + fac[i] + "=" + lvl[i]; 	
	xv.fac2lvl[ fac2 ] = lvl2;	
      }
    
    xv.name = xname;
    xv.type = var.type;
    xv.label = xlabel + ")";
    xv.base = var.name;
    xv.domain_name = var.domain_name;
    xv.domain_group = var.domain_group;

    // change to numeric encoding instead (but keeping all other info
    // for the data dictionary?)
    if ( options.numeric_strata_encoding )
      {
	int n = strata2number.size();
	if ( strata2number.find( xv.name ) == strata2number.end() )
	  strata2number[ xv.name ] = n+1;

	std::stringstream ss;
	ss << var.name << "." << strata2number[  xv.name ] ;
	xv.name = ss.str();
      }
    
    // store
    xvars.insert( xv );      

    return xv;
  }

  void check_variable_lengths();
  
  void check_variables_across_domains();

  // safe to return a const pointer to domain set elements, which are never
  // modified at the study data stage...
  const domain_t * domain( const std::string & name , const std::string & group ) const
  {
    std::set<domain_t>::const_iterator ii = domains.find( domain_t( name , group ) );
    if ( ii == domains.end() ) return NULL;
    return &(*ii);
  }

  void add( const indiv_t & indiv )
  {
    std::set<indiv_t>::iterator ii = indivs.find( indiv ) ;
    // create from this / or add to existing
    if ( ii == indivs.end() )
      {
	indivs.insert( indiv );
      }
    else
      {
	indiv_t temp = *ii;
	indivs.erase( ii );
	temp.add( indiv );
	indivs.insert( temp );
      }
  }
};

#endif
