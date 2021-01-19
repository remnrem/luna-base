

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

#include "merge.h"
#include "merge-helpers.h"
#include <fstream>
#include <cstdlib>
#include <wordexp.h>


// helpers


std::string expand( const std::string & f )
{
  wordexp_t exp_result;
  wordexp( f.c_str(), &exp_result, 0);
  std::string r = exp_result.we_wordv[0];
  wordfree(&exp_result);
  return r;
}



bool fileExists( const std::string & f )
{
  
  std::ifstream inp;
  
  inp.open(f.c_str(), std::ifstream::in);
  
  if( inp.fail() )
    {
      inp.clear(std::ios::failbit);
      inp.close();
      return false;
    }
  inp.close();
  return true;
  
}

bool file_extension( const std::string & f, const std::string & ext )
{
  int l = ext.size() + 1;  
  if ( f.size() < l ) return false;
  const std::string s = f.substr( f.size() - l );
  return iequals( s , "." + ext );
}

std::string remove_extension( const std::string & f, const std::string & ext )
{
  if ( ! file_extension( f, ext ) ) return f;
  return f.substr( 0 , f.size() - ( ext.size() + 1 ) ) ;  
}

bool iequals(const std::string& a, const std::string& b)
{
  unsigned int sz = a.size();
  if (b.size() != sz)
    return false;
  for (unsigned int i = 0; i < sz; ++i)
    if (tolower(a[i]) != tolower(b[i]))
      return false;
  return true;
}



std::string toupper( const std::string & s )
{
  std::string j = s;
  for (int i=0;i<j.size();i++) j[i] = std::toupper( s[i] );
  return j;
}


std::string sanitize( const std::string & s )
{
  std::string j = s;
  for (int i=0;i<j.size();i++)
    {
      if ( j[i] == '-' ) j[i] = '_';
      if ( j[i] == '+' ) j[i] = '_';
      if ( j[i] == ' ' ) j[i] = '_';
      if ( j[i] == '/' ) j[i] = '_';
      if ( j[i] == '\\' ) j[i] = '_';
      if ( j[i] == '*' ) j[i] = '_';
      if ( j[i] == '<' ) j[i] = '_';
      if ( j[i] == '>' ) j[i] = '_';
      if ( j[i] == '=' ) j[i] = '_';
      if ( j[i] == '&' ) j[i] = '_';
      if ( j[i] == '^' ) j[i] = '_';
      if ( j[i] == '!' ) j[i] = '_';
      if ( j[i] == '@' ) j[i] = '_';
      if ( j[i] == '#' ) j[i] = '_';
      if ( j[i] == '$' ) j[i] = '_';
      if ( j[i] == '%' ) j[i] = '_';
      if ( j[i] == '(' ) j[i] = '_';
      if ( j[i] == ')' ) j[i] = '_';     
    }
  return j;  
}



std::vector<std::string> parse(const std::string & item, const std::string & s , bool empty )
{  
  if ( s.size() == 1 ) return char_split( item , s[0] , empty ); 
  if ( s.size() == 2 ) return char_split( item , s[0] , s[1] , empty ); 
  if ( s.size() == 3 ) return char_split( item , s[0] , s[1] , s[2] , empty ); 
  halt("silly internal error in parse/char_split");
  std::vector<std::string> dummy;
  return dummy;
}  


std::vector<std::string> char_split( const std::string & s , const char c , bool empty )
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

std::vector<std::string> char_split( const std::string & s , const char c , const char c2 , bool empty )
{
  std::vector<std::string> strs;  
  if ( s.size() == 0 ) return strs;
  int p=0;

  for (int j=0; j<s.size(); j++)
    {	        
      if ( s[j] == c || s[j] == c2 ) 
	{ 	      
	  if ( j == p ) // empty slot?
	    {
	      if (empty) strs.push_back( "." );
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


std::vector<std::string> char_split( const std::string & s , const char c , const char c2 , const char c3 , bool empty )
{
  std::vector<std::string> strs;  
  if ( s.size() == 0 ) return strs;
  int p=0;

  for (int j=0; j<s.size(); j++)
    {	        
      if ( s[j] == c || s[j] == c2 || s[j] == c3 ) 
	{ 	      
	  if ( j == p ) // empty slot?
	    {
	      if (empty) strs.push_back( "." );
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


void halt( const std::string & msg )
{ 
  std::cerr << "\n*** error : " << msg << "\n";   
  std::exit(1);
}


std::istream& safe_getline(std::istream& is, std::string& t)
 {
   t.clear();
   
   // The characters in the stream are read one-by-one using a std::streambuf.
   // That is faster than reading them one-by-one using the std::istream.
   // Code that uses streambuf this way must be guarded by a sentry object.
   // The sentry object performs various tasks,
   // such as thread synchronization and updating the stream state.
   
   std::istream::sentry se(is, true);
   std::streambuf* sb = is.rdbuf();
   
   for ( ; ; ) 
     {
       
       int c = sb->sbumpc();
       
       switch (c) 
	 {
	 case '\n':
	   return is;
	   
	 case '\r':
	   if (sb->sgetc() == '\n')
	     sb->sbumpc();
	   return is;
	   
	   // replace w/ macro EOF to compile	  
	   // 	case std::streambuf::traits_type::eof() :
	 case EOF :
	   // Also handle the case when the last line has no line ending
	   if(t.empty())
	     is.setstate(std::ios::eofbit);
	   return is;
	   
	 default:
	   t += (char)c;
	 }
     }
 }

bool imatch(const std::string& a, const std::string& b , unsigned int min )
{
  // only compare up to length of the shortest;
  // if min specified (default 0) then compare at least min chars

  unsigned int sz = a.size() < b.size() ? a.size() : b.size() ;
  if ( min != 0 ) sz = min;
  if ( a.size() < min || b.size() < min ) return false;
  
  for (unsigned int i = 0; i < sz; ++i)
    if (tolower(a[i]) != tolower(b[i]))
      return false;
  return true;
}

bool yesno( const std::string & s )
{
  // 0 no NO n N F f false FALSE 
  // versus all else 
  if ( s.size() == 0 ) return false; // empty == NO
  if ( s[0] == '0' || s[0] == 'n' || s[0] == 'N' || s[0] == 'f' || s[0] == 'F' ) return false;
  return true;
}

std::string search_replace( const std::string & s , char a , char b )
{
  std::string j = s;
  for (int i=0;i<j.size();i++) if ( j[i] == a ) j[i] = b;
  return j;
}


bool str2int(const std::string & s , int * i)
{
  return from_string<int>(*i,s,std::dec);
}

bool str2dbl(const std::string & s , double * d)
{
  return from_string<double>(*d,s,std::dec);
}
