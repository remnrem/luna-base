
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
#include <map>
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


bool is_missing( const std::string & s )
{
  return s == "NA" || s == "." || s == "";
}

void usage()
{
  std::cerr << "usage: fixrows [--keep-all-missing|-k] KEY_COL [KEY_COL ...]\n";
}


struct row_t { 
  row_t( const std::vector<std::string> & h ,
	 const std::vector<std::string> & d ,
	 const std::set<std::string> & rows ,
	 std::map<std::string,std::string> * vars )
  {

    // populate either as a row-identifier, or a variable
    for (int i=0;i<h.size();i++)
      if ( rows.find( h[i] ) != rows.end() )
	faclvl[ h[i] ] = d[i] ;
      else
	(*vars)[ h[i] ] = d[i];
  } 
  
  std::map<std::string,std::string> faclvl;
  
  bool operator<( const row_t & rhs ) const { 
    if ( faclvl.size() < rhs.faclvl.size() ) return true;
    if ( faclvl.size() > rhs.faclvl.size() ) return false;
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


void addin( std::map<std::string,std::string> * m , const std::map<std::string,std::string> & n )
{
  std::map<std::string,std::string>::const_iterator nn = n.begin();
  while ( nn != n.end() )
    {
      // add if not present, no matter what the new value is
      if ( m->find( nn->first ) == m->end() )
	{
	  (*m)[ nn->first ] = nn->second;
	}
      else // otherwise, only add non-missing values (but check for conflict)
	{

	  // if N missing, then ignore
	  if ( is_missing( nn->second ) ) { ++nn; continue; }

	  // if M missing, then add
	  if ( is_missing( (*m)[ nn->first ] ) )
	    (*m)[ nn->first ] = nn->second;	    
	  else
	    {
	      // otherwise, just check that we did not have different values
	      
	      if ( nn->second != (*m)[ nn->first ] )
		{
		  std::cerr << " problem, found non-identical non-missing values\n";
		  std::cout << nn->first << " " << nn->second << " " << (*m)[ nn->first ] << "\n";
		  std::exit(1);
		}
	    }	  
	}     
      ++nn;
    }
}


int main(int argc , char ** argv )
{
  
  // ./fixrows ID CH < old.txt > new.txt

  // opts = factors to fix / make unique
  //        merge values, replacing 'NA' with the observed value
  //        write out again

  // datastore:
  // id --> variable --> value 
  std::map<row_t,std::map<std::string,std::string> > store;

  // default: drop non-key columns that are entirely missing post-merge
  bool drop_all_missing = true;
  static struct option long_options[] =
    {
      {"keep-all-missing", no_argument, 0, 'k'},
      {0, 0, 0, 0}
    };

  while ( true )
    {
      int option_index = 0;
      int c = getopt_long( argc, argv, "k", long_options, &option_index );
      if ( c == -1 ) break;
      if ( c == 'k' )
        drop_all_missing = false;
      else
        {
          std::cerr << "fixrows: unrecognized option\n";
          usage();
          std::exit(1);
        }
    }
  
  // which factors to uniqify
  std::set<std::string> facs;
  for (int i=optind;i<argc;i++)
    {
      std::string a = argv[i];
      facs.insert( a ); 
    }
  
  // get all input from STDIN, read into memory 
  // if col 1 == "ID" this is a header

  std::vector<std::string> headers;
  std::vector<std::string> data;

  //
  // Input
  //

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
	  
	  if ( data.size() != col ) 
	    {
	      std::cerr << "fixrows: read bad line...\n";
	      std::exit(0);
	    }
	      
	  // add to store: get row key
	  // get variable names (based on header)
	  std::map<std::string,std::string> vars;
	  row_t row( headers , data , facs , &vars );

	  // merge, respecting NA's
	  addin( &( store[ row ] ) , vars );
	  
	}
    }
  
  
  //
  // output header
  //

  std::vector<int> keep_cols;
  keep_cols.reserve( col );

  for (int i=0; i<col; i++)
    {
      const std::string & h = headers[i];
      bool keep = true;

      if ( drop_all_missing && facs.find( h ) == facs.end() )
        {
          keep = false;
          std::map<row_t,std::map<std::string,std::string> >::const_iterator tt = store.begin();
          while ( tt != store.end() )
            {
              std::map<std::string,std::string>::const_iterator vv = tt->second.find( h );
              const std::string value = vv == tt->second.end() ? "" : vv->second;
              if ( ! is_missing( value ) )
                {
                  keep = true;
                  break;
                }
              ++tt;
            }
        }

      if ( keep ) keep_cols.push_back( i );
    }

  for (int j=0; j<keep_cols.size(); j++)
    {
      const int i = keep_cols[j];
      std::cout << ( j != 0 ? "\t" : "" ) << headers[i] ;
    }
  std::cout << "\n";

  //
  // output data rows
  //

  std::map<row_t,std::map<std::string,std::string> >::const_iterator ss = store.begin();
  while ( ss != store.end() )
    {
      const std::map<std::string,std::string> & vars = ss->second;
      
      for (int j=0; j<keep_cols.size(); j++)
	{
	  const int i = keep_cols[j];
	  const std::string & v = headers[i];

	  if ( j != 0 ) std::cout << "\t";
	  
	  if ( vars.find( v ) != vars.end() ) 
	    std::cout << vars.find( v )->second;
	  else
	    {	      
	      std::map<std::string,std::string>::const_iterator vv = ss->first.faclvl.find( v );
	      if ( vv == ss->first.faclvl.end() )
		{
		  std::cerr << "fixrows: internal error, this should not happen...\n";
		  std::exit(1);
		}
	      std::cout << vv->second;
	    }
	}      
      std::cout << "\n";
      ++ss;      
    }
  exit(0);
}
