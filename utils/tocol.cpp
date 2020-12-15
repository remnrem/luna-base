
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
#include <cstdlib>

#include <vector>
#include <cstring>
#include <sstream>

#include <unistd.h>
#include <getopt.h>
#include <cstdlib>

#include <map>
#include <set>

const int MAXBUF = 50000;

std::vector<std::string> char_split( const std::string & s , const char c , bool empty = true )
{

  std::vector<std::string> strs;  
  if ( s.size() == 0 ) return strs;
  int p=0;

  for (int j=0; j<s.size(); j++)
    {	        
      if ( s[j] == c ) 
	{ 	      
	  if ( j == p ) // empty slot?
	    {
	      if ( empty ) strs.push_back( "." );
	      ++p;
	    }
	  else
	    {
	      strs.push_back(s.substr(p,j-p)); 
	      p=j+1; 
	    }
	}	  
    }
  
  if ( empty && p == s.size() ) 
    strs.push_back( "." );
  else if ( p < s.size() )
    strs.push_back( s.substr(p) );
  
  return strs;
}


bool tokenize( std::vector<std::string> & data , const int n )
{  
  char line[MAXBUF];
  std::cin.getline( line, MAXBUF, '\n' );  
  std::string sline = line;  
  data = char_split( sline , '\t' );  
  return data.size() == n; 
}


struct row_t { 
  row_t( const std::vector<std::string> & h ,
	 const std::vector<std::string> & d ,
	 const std::set<std::string> & rows ,
	 const std::set<std::string> & cols ,
	 std::map<std::string,std::string> * vars )
  {
    // ID always first column
    id = d[0];
    // rows
    for (int i=1;i<h.size();i++)
      if ( rows.find( h[i] ) != rows.end() )
	faclvl[ h[i] ] = d[i] ;
    
    // col str
    std::string cstr;
    for (int i=1;i<h.size();i++)
      if ( cols.find( h[i] ) != cols.end() )
	cstr += "." + h[i] + "_" + d[i];
    
    // populate var names, minus factors 
    for (int i=1;i<h.size();i++)
      if ( rows.find( h[i] ) == rows.end() && cols.find( h[i] ) == cols.end() )
	(*vars)[ h[i] + cstr ] = d[i];

  } 
	 
  std::string id;
  std::map<std::string,std::string> faclvl;
  bool operator<( const row_t & rhs ) const { 
    if ( id < rhs.id ) return true;
    if ( id > rhs.id ) return false;
    if ( faclvl.size() < rhs.faclvl.size() ) return true;
    if ( faclvl.size() > rhs.faclvl.size() ) return true;
    std::map<std::string,std::string>::const_iterator ii = faclvl.begin();
    std::map<std::string,std::string>::const_iterator jj = rhs.faclvl.begin();
    while ( ii != faclvl.end() )
      {
	if ( ii->first < jj->first ) return true;
	if ( ii->first > jj->first ) return false;
	if ( ii->second < jj->second ) return true;
	if ( ii->second > jj->second ) return false;
	++ii;
	++jj;
      }
    return false;
  }
};



int main(int argc , char ** argv )
{
  
  // tocol fac1 fac2 < input > output 
  
  // take long format across multiple files; 
  // all files expected to 1) have a header, 2) have the same header
  // make column format (i.e. 

  // id --> variable --> value 
  std::map<row_t,std::map<std::string,std::string> > store;
  
  // which factors (columns) to make wide?
  bool torow = true;
  std::set<std::string> rows, cols;
  for (int i=1;i<argc;i++)
    {
      std::string a = argv[i];
      if ( a == "/" ) torow = false;
      if ( torow ) rows.insert( a ); 
      else cols.insert( a );
    }

  // get all input from STDIN, read into memory 
  // if col 1 == "ID" this is a header
  std::vector<std::string> headers;
  std::vector<std::string> data;

  int col = 0;

  while ( ! std::cin.eof() ) 
    {

      // get line
      char line[MAXBUF];
      std::cin.getline( line, MAXBUF, '\n' );  
      std::string sline = line;  

      // first header?
      if ( col == 0 ) 
	{
	  headers = char_split( sline , '\t' );  
	  col = headers.size();
	  data.resize( col );
	}
      else
	{
	  data = char_split( sline , '\t' );
	  if ( data.size() == 0 ) continue;
	  
	  // a new header?
	  if ( data[0] == "ID" ) 
	    {
	      headers = data;
	      col = headers.size();
	      data.resize( col );
	    }
	  else // parse as a data line
	    {
	      if ( data.size() != col ) 
		{
		  std::cerr << "read bad line...\n";
		  std::exit(0);
		}
	      
	      // add to store: get row key
	      // get variable names (based on header plus cols
	      std::map<std::string,std::string> vars;
	      row_t row( headers , data , rows , cols , &vars );	      
	      
	      store[ row ] =  vars ;
	    
	    }
	}
    }

  // now output
  
  std::cout << "read " << store.size() << " rows\n";

  exit(0);
}

