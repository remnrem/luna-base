
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

enum type_t { FACTOR , TEXT , INT , FLOAT , YESNO , DATE , TIME };
bool type_check( const std::string & value , type_t );


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
    missing = "";
  }
  
  domain_t( const std::string & filename )
  {
    read( filename );
  }
  
  int read( const std::string & filename );  

  std::string name; // e.g. 'eeg'
  std::string group; // e.g. 'spindles'

  std::map<std::string,var_t> variables;
  
  // domain-specific missing data symbol?
  std::string missing;
  
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
    return group < rhs.name ; 
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
	  halt( "multiple values for " + id + " " + ii->first.name );
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
  
  std::set<domain_t> domains;

  std::set<var_t> xvars; // expanded vars, e.g. DENS.F_11_CH_C3_SS_N2

  std::map<std::string,std::string> faclabels; // track labels for factors (for dict output)
  
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
    // already done expansion?
    //
    
    if ( xvars.find( var ) != xvars.end() )
      return var;

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
	xname += (i > 0 ? "_" : "." ) + fac2 + "_" + lvl2;
	xlabel += (i > 0 ? ", " : "" ) + fac[i] + "=" + lvl[i]; 	
	xv.fac2lvl[ fac2 ] = lvl2;	
      }
    
    xv.name = xname;
    xv.type = var.type;
    xv.label = xlabel + ")";
    xv.base = var.name;
    xv.domain_name = var.domain_name;
    xv.domain_group = var.domain_group;

    // store
    xvars.insert( xv );      

    return xv;
  }
  
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
