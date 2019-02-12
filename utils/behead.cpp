#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <vector>
#include <cstring>
#include <sstream>

#include <unistd.h>
#include <getopt.h>
#include <cstdlib>

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

int main(int argc , char ** argv )
{
  
  bool print_nums = false;
  bool tab_sep = false;
  
  for (int i=1;i<argc;i++)
    {
      
      if      ( strcmp( argv[i] , "-n" ) == 0 ) print_nums = true;
      else if ( strcmp( argv[i] , "-t" ) == 0 ) tab_sep = true;
      else if ( strcmp( argv[i] , "-nt" ) == 0 ) { tab_sep = true; print_nums = true;}
      else if ( strcmp( argv[i] , "-tn" ) == 0 ) { tab_sep = true; print_nums = true;}
      else { std::cerr << "did not recognize option " << argv[i] << "\n"; std::exit(1); } 
    }      


  int col = 0;
  int row = 0;

  std::vector<std::string> headers;
  std::vector<std::string> data;

  while ( ! std::cin.eof() ) 
    {
      // parse a line?
      if ( col ) 
	{
	  
	  if ( ! tokenize( data , col ) ) 
	    {
	      if ( data.size() > 0 ) std::cout << "read bad line...\n";
	      continue;
	    }
	  
	  ++row;
	  
	  for (int i=0; i<col; i++ ) 
	    {

	      if ( print_nums ) 
		{
		  if ( tab_sep  )
		    std::cout << row << "\t" << i+1 << "\t";
		  else 
		    std::cout << std::setw(6) << row << std::setw(6) << i+1 ;
		}
	      
	      if ( tab_sep )
		std::cout << headers[i] << "\t" << data[i] << "\n";
	      else
		std::cout << std::right << std::setw(25) << headers[i] 
			  << "   "
			  << std::left << std::setw(20)  << data[i] << "\n";
	    }
	  
	  // space between each individual
	  std::cout << "\n";
	}
      else // this is first row -- read as headers
	{
	  char line[MAXBUF];
	  std::cin.getline( line, MAXBUF, '\n' );
	  std::string buf;
	  std::string sline = line;
	  headers = char_split( sline , '\t' );
	  col = headers.size();
	  data.resize( col );
	}
	
    }
  exit(0);
}

