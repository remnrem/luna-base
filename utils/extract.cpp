
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
#include <map>
#include <set>
#include <string>
#include <vector>
#include "helper/helper.h"

struct row_t
{


  row_t( const std::string & id , 
	 const std::string & tag )
    : id(id) , tag(tag) { } 

  row_t( const std::string & id , 
	 const std::string & tag , 
	 const std::vector<std::string> & lnames ,
	 const std::vector<std::string> & lvalues ) : 
    id(id) , tag(tag) , lnames(lnames) , lvalues(lvalues)
  {    
  }

  std::string id;
  std::string tag;
  std::vector<std::string> lnames;
  std::vector<std::string> lvalues;
  
  bool operator<( const row_t & rhs ) const
  {
    if ( id < rhs.id ) return true;
    if ( id > rhs.id ) return false;
    if ( tag < rhs.tag ) return true;
    if ( tag > rhs.tag ) return false;
    if ( lvalues.size() < rhs.lvalues.size() ) return true;
    if ( lvalues.size() > rhs.lvalues.size() ) return false;
    for (int j=0;j<lvalues.size();j++) 
      {
	if ( lvalues[j] < rhs.lvalues[j] ) return true;
	if ( lvalues[j] > rhs.lvalues[j] ) return false;
      }
    return false;
  }
};


struct data_t
{  
  static std::set<std::string> allvar;
  data_t() { d.clear(); }  
  void add( const std::string & var , const std::string & val ) 
  { 
    // check this hasn't been seen before
    if ( d.find( var ) != d.end() )
      {
	if ( d[var] == val )
	  std::cerr << "duplicated row, same value\n";
	else
	  Helper::halt("duplicated row, different values");
      }
    allvar.insert(var);
    d[var] = val; 
  }  
  std::map<std::string,std::string> d;
};

std::set<std::string> data_t::allvar;

int main( int argc , char ** argv )
{
  
  //
  // Format
  //
  
  // extract -i=@id.txt -v=@v.txt -klevel=x -dlevel=y < input.long > output.dat
  
  // -i  extract only these individuals (ID)
  // -v  extract only these variables
  // -t  keep onlt these tags
  
  // -k  keep only these levels  
  // -d  keep only these values
  // -c  make cols by this level
  
  // input  : ID  tag var {level1} {level2} {value} 
  
  //  levels have format level=value
  //   we can either collapse by levels
  //     e.g.   ID  V1_level1_x  V1_level1_y V2_level1_z
  //   or we can keep level as a row-level item
  //     e.g.   ID  level1   V1
  //            id1 x        .
  //            id1 y        .
  //            id1 z        .
  //            id2 x        .
  //   and/or we can drop levels 
  
  // output : row (default) by ID,tag 
  //        : row (alt)  ID, level1
  //        : row (alt)  ID, level1
  //        : row (alt)  ID, level-x,level-y,...
  
  // options: can collapse tag-value onto the same level
  
  std::set<std::string> xindiv;
  std::set<std::string> xvar;
  std::set<std::string> xtag;
  std::set<std::string> xlevel;
  std::set<std::string> dlevel;
  std::set<std::string> clevel;
  
  // process command line
  if ( argc > 1 )
    {
      for (int i=1; i < argc; i++) 
	{
	  std::string s = argv[i];
	  std::cout << "s [" << s << "]\n";
	  if ( s.size() < 4 || s[0] != '-' || s[2] != '=' ) 
	    Helper::halt("expecting cmd args: -i, -v, -t, -c, -d or -k");
	  
	  char t = argv[i][1];
	  if ( t != 'i' && t != 'v' && t != 'd' && t != 'k' && t != 't' )
	    Helper::halt("expecting cmd args: -i, -v, -t, -c, -d or -k");
	  
	  s = s.substr( 3 );
	  
	  // reading from a file?  @includes
	  bool from_file = s[0] == '@';
	  if ( from_file ) s = s.substr(1);
	  
	  if ( from_file )
	    {
	    }
	  else
	    {
	      std::vector<std::string> tok = Helper::parse( s , "," );
	      for (int j=0;j<tok.size();j++) 
		{
		  if      ( t == 'i' ) xindiv.insert( tok[j] );
		  else if ( t == 'v' ) xvar.insert( tok[j] );
		  else if ( t == 't' ) xtag.insert( tok[j] );
		  else if ( t == 'k' ) xlevel.insert( tok[j] );
		  else if ( t == 'd' ) dlevel.insert( tok[j] );
		}
	    }
	  
	}
    }
  
  bool filter_indiv  = xindiv.size() > 0 ;
  bool filter_tag    = xtag.size() > 0 ;
  bool filter_var    = xvar.size() > 0 ;
  bool filter_klevel = xlevel.size() > 0 ;
  bool filter_dlevel = dlevel.size() > 0 ;
  bool collapse_level = clevel.size() > 0 ;
  
  if ( filter_klevel && filter_dlevel ) Helper::halt( "cannot specifify both -d and -k" );
  
  // extract -i=@id.txt -v=@v.txt -klevel=x -dlevel=y < input.long > output.dat

  // 
  // Store all data
  //
  
  std::map<row_t,data_t> data;
  

  //
  // Read all input 
  //
  
  
  
  while ( ! std::cin.eof() )
    {
      
      std::string line;
      std::getline( std::cin , line);
      
      if ( std::cin.eof() ) break;
      if ( line == "" ) continue;
      
      // ID tag var (l1=s1} {l2=s2} value
      std::vector<std::string> tok = Helper::parse( line , "\t" );
      if ( tok.size() == 0 ) continue;
      if ( tok.size() < 4 ) Helper::halt( "bad line: " + line );
      int nlevels = tok.size() - 4;
      
      std::string id = tok[0];
      std::string tag = tok[1];
      std::string var = tok[2];

      if ( filter_indiv && xindiv.find( id ) == xindiv.end() ) 
	continue;

      if ( filter_tag && xtag.find( tag ) == xtag.end() ) 
	continue;
      
      if ( filter_var && xvar.find( var ) == xvar.end() ) 
	continue;
      
      // levels?
      std::vector<std::string> col_levels, row_level_names, row_level_values;

      for ( int l = 0 ; l < nlevels ; l++ )
	{
	  std::string & lvl = tok[ l + 3 ];
	  
	  int lsplit = lvl.find( "=" );
	  if ( lsplit == std::string::npos ) 
	    Helper::halt( "bad level format, expecting level=value" );
	  
	  std::string level_name = lvl.substr( 0 , lsplit );
	  std::string level_value = lvl.substr( lsplit + 1 );
	  
	  // drop this level?
	  if ( filter_dlevel && dlevel.find( lvl ) != dlevel.end() )
	    continue;
	  if ( filter_klevel && xlevel.find( lvl ) == xlevel.end() )
	    continue;

	  // this level as a row or col?
	  
	  bool ascol = collapse_level && clevel.find( level_name ) != clevel.end() ;
	  
	  if ( ascol ) col_levels.push_back( level_name + "." + level_value );
	  else 
	    {
	      row_level_names.push_back( level_name );
	      row_level_values.push_back( level_value );
	    }
	}
      
      
      //
      // the actual value
      //
      
      std::string val = tok[ tok.size() - 1 ];
      
      
      //
      // if we are collapsing levels, we are adding multiple (i.e. creating new) 
      // variables; otherwise we are adding just a single  variable
      //
      
      row_t row( id , tag , row_level_names , row_level_values ); 
      
      // collapsing levels? 
      
      if ( collapse_level && col_levels.size() > 0 ) 
	{
	  for (int j = 0 ; j < col_levels.size() ; j++ )
	    {
	      std::string new_var = var + "." + col_levels[j];
	      data[ row ].add( new_var , val );
	    }
	}
      else
	{	  
	  data[ row ].add( var , val );
	}
    }

  //
  // Done reading all input; find the total number of columns that we will need to output
  //
  
  std::cerr << "will output " << data.size() << " rows and " << data_t::allvar.size() << " cols\n";

  //
  // Now output everything
  //

  // header
  
  std::cout << "ID" 
	    << "\tTAG";
  
  // levels

  std::set<std::string>::const_iterator cc = data_t::allvar.begin();
  while ( cc != data_t::allvar.end() )
    {
      std::cout << "\t" << *cc ;
      ++cc;
    }
  std::cout << "\n";
  
  // each row/individual
  
  std::map<row_t,data_t>::const_iterator rr = data.begin();
  while ( rr != data.end() )
    {
      
      const row_t & r = rr->first;
      
      std::cout << r.id 
		<< "\t" << r.tag;
      
      const data_t & d = rr->second;
      
      std::set<std::string>::const_iterator cc = data_t::allvar.begin();
      while ( cc != data_t::allvar.end() )
	{
	  std::map<std::string,std::string>::const_iterator ff = d.d.find(*cc);
	  if ( ff != d.d.end() ) std::cout << "\t" << ff->second;
	  else std::cout << "\t.";
	  ++cc;
	}
      std::cout << "\n";
      
      ++rr;

    }

  std::exit(0);
}
