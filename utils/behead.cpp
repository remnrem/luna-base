#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <vector>
#include <cstring>
#include <string>

std::vector<std::string> char_split( const std::string & s , const char c , bool empty = true )
{

  std::vector<std::string> strs;  
  if ( s.size() == 0 ) return strs;
  std::string::size_type p = 0;

  for ( std::string::size_type j = 0 ; j < s.size() ; ++j )
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


bool tokenize( std::vector<std::string> & data , const std::string & line , const int n )
{  
  data = char_split( line , '\t' );  
  return static_cast<int>( data.size() ) == n; 
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
  std::string line;

  if ( ! std::getline( std::cin , line ) ) return 0;

  headers = char_split( line , '\t' );
  col = headers.size();
  data.resize( col );

  while ( std::getline( std::cin , line ) ) 
    {
      
      if ( ! tokenize( data , line , col ) ) 
	{
	  if ( data.size() > 0 ) std::cerr << "read bad line...\n";
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
  exit(0);
}
