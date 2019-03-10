
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

#include "helper.h"
#include "logger.h"

#include "defs/defs.h"
#include "intervals/intervals.h"

#include <cmath>
#include <cstdio>

#include <iostream>
#include <iomanip>
#include <fstream>

#ifndef WINDOWS
#include <wordexp.h>
#endif

extern logger_t logger;

std::string Helper::toupper( const std::string & s )
{
  std::string j = s;
  for (int i=0;i<j.size();i++) j[i] = std::toupper( s[i] );
  return j;
}


std::string Helper::sanitize( const std::string & s )
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

std::string Helper::search_replace( const std::string & s , char a , char b )
{
  std::string j = s;
  for (int i=0;i<j.size();i++) if ( j[i] == a ) j[i] = b;
  return j;
}

std::string Helper::expand( const std::string & f )
{
#ifdef WINDOWS
  return f;
#else
  wordexp_t exp_result;
  wordexp( f.c_str(), &exp_result, 0);
  std::string r = exp_result.we_wordv[0];
  wordfree(&exp_result);
  return r;
#endif
}


bool Helper::is_folder( const std::string & f ) { if ( f.size() == 0 ) return false; return f[f.size()-1]== globals::folder_delimiter; } 


bool Helper::file_extension( const std::string & f, const std::string & ext )
{
  int l = ext.size() + 1;  
  if ( f.size() < l ) return false;
  const std::string s = f.substr( f.size() - l );
  return Helper::iequals( s , "." + ext );
}

void Helper::halt( const std::string & msg )
{
  
  // some other code handles the exit, e.g. if running under luna-web?
  if ( globals::bail_function != NULL ) 
    globals::bail_function( msg );
  
  // switch logger off , i.e. as we don't want close-out msg
  logger.off();

  // generic bail function (not using logger)
  std::cerr << "error : " << msg << "\n";   
  std::exit(1);
}

void Helper::warn( const std::string & msg )
{
  logger.warning( msg );
}

void Helper::debug( const std::string & msg )
{
  std::cerr << "debug : " << msg << "\n";
}


bool Helper::realnum(double d)
{
  double zero = 0;
  if (d != d || d == 1/zero || d == -1/zero)
    return false;
  else
    return true;
}


bool Helper::similar( double a, double b , double EPS )
{
  return fabs( a - b ) < EPS ;
}


std::string Helper::int2str(int n)
{
  std::ostringstream s2( std::stringstream::out );
  s2 << n;
  return s2.str();
}

std::string Helper::int2str(long n)
{
  std::ostringstream s2( std::stringstream::out );
  s2 << n;
  return s2.str();
}

std::string Helper::int2str(uint64_t n)
{
  std::ostringstream s2( std::stringstream::out );
  s2 << n;
  return s2.str();
}

std::string Helper::dbl2str(double n)
{
  std::ostringstream s2( std::stringstream::out );
  s2 << n;
  return s2.str();
}

std::string Helper::dbl2str(double n, int dp )
{
  std::ostringstream ss( std::stringstream::out );
  ss << std::fixed
     << std::setprecision( dp );
  ss << n;
  return ss.str();
}

std::string Helper::dbl2str_fixed(double n, int ch )
{

  //
  // first, check if a standard conversion works
  //

  std::ostringstream s2( std::stringstream::out );
  s2 << n;
  std::string retstr1 = s2.str();
  if ( retstr1.size() <= ch ) return retstr1;
  
  // if not, need to try to use a better encoding
  // that allows total of 'ch' characters

  double a = n > 0 ? n : -n;
  int ch2 = n < 0 ? ch-1 : ch; // need to allow for negative sign

  // with 8 characters, max value is 10e8 - 1 = 99,999,999 (+ve)
  //                                 10e7 - 1 = -9,999,999 (-ive)

  if ( a >= pow( 10 , ch2 ) ) 
    {
      std::cerr << "trying to print " << n << " in " << ch << " characters...\n";
      Helper::halt( "EDF silliness: need to rescale signal channel so min/max can be represented in 8 chars...");
    }

  // use fixed format; take as many DP as possible. 
  
  std::ostringstream s3( std::stringstream::out );
  s3 << std::fixed
     << std::setprecision( ch );
  s3 << n;
  std::string retstr = s3.str();
  return retstr.substr(0,ch);

}

bool Helper::str2dbl(const std::string & s , double * d)
{
  return from_string<double>(*d,s,std::dec);
}

bool Helper::str2int(const std::string & s , int * i)
{
  return from_string<int>(*i,s,std::dec);
}

bool Helper::str2int64(const std::string & s , uint64_t * i)
{
  return from_string<uint64_t>(*i,s,std::dec);
}


// std::string Helper::stringize( const std::set<std::string> & d , const std::string & delim )
// {
//   std::stringstream ss;
//   std::set<std::string>::const_iterator dd = d.begin();
//   while ( dd != d.end() )
//     {
//       if ( dd != d.begin() ) ss << delim;
//       ss << *dd;
//       ++dd;
//     }
//   return ss.str();
// }

std::string Helper::trim( const std::string & s , const char c , const char d )
{
  int first_nonspace = 0;
  int last_nonspace = s.size()-1;
  for (int i=0;i<s.size();i++) if ( s[i] == c || s[i] == d ) ++first_nonspace; else break;
  for (int i=s.size()-1;i!=0;i--) if ( s[i] == c || s[i] == d ) --last_nonspace; else break;  
  if ( first_nonspace > last_nonspace ) return "";
  return s.substr( first_nonspace , last_nonspace - first_nonspace  + 1 );
}

std::string Helper::format( const std::string & s , int indent , int width , bool no_initial_indent )
{
  std::string r;  
  int p = 0;  
  while ( p < s.size() )
    {      
      if ( p > 0 || ! no_initial_indent )
	r += std::string( indent , ' ' );      
      int x = p + width - indent;
      if ( x >= s.size() )
	{
	  r += s.substr( p );
	  break;
	}
      else
	{
	  // find next space
	  while (x>p && s[x] != ' ' && s[x] != '-' ) { --x; }
	  r += s.substr( p , x - p ) + "\n";
	  p = x+1;
	}
    }
  return r;  
}


std::vector<std::string> Helper::parse(const std::string & item, const std::string & s , bool empty )
{  
  if ( s.size() == 1 ) return Helper::char_split( item , s[0] , empty ); 
  if ( s.size() == 2 ) return Helper::char_split( item , s[0] , s[1] , empty ); 
  if ( s.size() == 3 ) return Helper::char_split( item , s[0] , s[1] , s[2] , empty ); 
  Helper::halt("silly internal error in parse/char_split");
  std::vector<std::string> dummy;
  return dummy;
}  

std::vector<std::string> Helper::quoted_parse(const std::string & item , const std::string & s , const char q , const char q2, bool empty )
{
  if ( s.size() == 1 ) return Helper::quoted_char_split( item , s[0] , q, q2, empty ); 
  if ( s.size() == 2 ) return Helper::quoted_char_split( item , s[0] , s[1] , q, q2, empty ); 
  if ( s.size() == 3 ) return Helper::quoted_char_split( item , s[0] , s[1] , s[2] , q, q2, empty ); 
  Helper::halt("silly internal error in parse/char_split");
  std::vector<std::string> dummy;
  return dummy;
}


std::vector<std::string> Helper::char_split( const std::string & s , const char c , bool empty )
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

std::vector<std::string> Helper::char_split( const std::string & s , const char c , const char c2 , bool empty )
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


std::vector<std::string> Helper::char_split( const std::string & s , const char c , const char c2 , const char c3 , bool empty )
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


std::vector<std::string> Helper::quoted_char_split( const std::string & s , const char c , const char q , const char q2, bool empty )
{

  std::vector<std::string> strs;  
  if ( s.size() == 0 ) return strs;
  int p=0;
  
  bool in_quote = false;
  
  for (int j=0; j<s.size(); j++)
    {	        

      if ( s[j] == '"' || s[j] == q || s[j] == q2 ) in_quote = ! in_quote;
      
      if ( (!in_quote) && s[j] == c ) 
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


std::vector<std::string> Helper::quoted_char_split( const std::string & s , const char c , const char c2 , const char q , const char q2 , bool empty )
{

  std::vector<std::string> strs;  
  if ( s.size() == 0 ) return strs;
  int p=0;

  bool in_quote = false;

  for (int j=0; j<s.size(); j++)
    {	        

      if ( s[j] == '"' || s[j] == q || s[j] == q2 ) in_quote = ! in_quote;

      if ( (!in_quote) && ( s[j] == c || s[j] == c2 ) ) 
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


std::vector<std::string> Helper::quoted_char_split( const std::string & s , const char c , const char c2 , const char c3 , const char q , const char q2, bool empty )
{

  std::vector<std::string> strs;  
  if ( s.size() == 0 ) return strs;
  int p=0;

  bool in_quote = false;

  for (int j=0; j<s.size(); j++)
    {	        

      if ( s[j] == '"' || s[j] == q || s[j] == q2 ) in_quote = ! in_quote;

      if ( (!in_quote) && ( s[j] == c || s[j] == c2 || s[j] == c3 ) ) 
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


std::string Helper::brief( const std::string & s , int l )
{
  if ( s.size() < l ) return s;
  return s.substr(0,l-3) + "...";
}



bool Helper::deleteFile( const std::string & f )
{
  if ( ! fileExists( f ) ) return false; 
  if ( remove( f.c_str() )  != 0 ) Helper::halt( "problem clearing database " + f );
  return true;
}

bool Helper::fileExists( const std::string & f )
{
  
  std::ifstream inp;
  
  inp.open(f.c_str(), std::ifstream::in);

  if(inp.fail())
    {
      inp.clear(std::ios::failbit);
      inp.close();
      return false;
    }
  inp.close();
  return true;

}


bool Helper::iequals(const std::string& a, const std::string& b)
{
  unsigned int sz = a.size();
  if (b.size() != sz)
    return false;
  for (unsigned int i = 0; i < sz; ++i)
    if (tolower(a[i]) != tolower(b[i]))
      return false;
  return true;
}


bool Helper::imatch(const std::string& a, const std::string& b)
{
  unsigned int sz = a.size() < b.size() ? a.size() : b.size() ;
  for (unsigned int i = 0; i < sz; ++i)
    if (tolower(a[i]) != tolower(b[i]))
      return false;
  return true;
}

bool Helper::yesno( const std::string & s )
{
  // 0 no NO n N F f false FALSE 
  // versus all else 
  if ( s.size() == 0 ) return false; // empty == NO
  if ( s[0] == '0' || s[0] == 'n' || s[0] == 'N' || s[0] == 'f' || s[0] == 'F' ) return false;
  return true;
}

clocktime_t::clocktime_t( const std::string & t )
{  
  valid = Helper::timestring( t , &h, &m, &s );
  if ( h < 0 || m < 0 || s < 0 ) valid = false;
  if ( h > 23 || m > 59 || s > 59 ) valid = false;  
};

bool clocktime_t::midpoint( const clocktime_t & t1 , const clocktime_t & t2 )
{

  if ( ! ( t1.valid && t2.valid ) )
    {
      valid = false;
      return false;
    }
 
  // copy first time
  h=t1.h; m=t1.m; s=t1.s;
  
  // time difference (from t1 to t2)
  // we assume that t1 is earlier, so +ve means 
  double diff = difference( t1 , t2 );
  
  // get midpoint
  diff /= 2.0;
  
  // advance this time by the required amount
  advance( diff );

  return true;
}

void clocktime_t::advance( uint64_t tp )
{
  // convert to hours
  double sec = tp / globals::tp_1sec;
  advance( sec / 3600.0 );
}



std::string Helper::timestring( const std::string & st , const interval_t & i )
{  
  int h0 = 0, m0 = 0, s0 = 0;  
  if ( ! Helper::timestring( st, &h0, &m0, &s0 ) ) return ".";  
  int h1 = h0, m1 = m0, s1 = s0;
  int h2 = h0, m2 = m0, s2 = s0;
  Helper::add_clocktime( &h1, &m1, &s1 , i.start );
  Helper::add_clocktime( &h2, &m2, &s2 , i.stop ); 
  
  std::stringstream str;
  str.precision(0);
  str << std::fixed;
  str << Helper::timestring(h1,m1,s1) << " - " << Helper::timestring(h2,m2,s2);
  return str.str(); 
}


std::string Helper::timestring( uint64_t a )
{

  // a is tp units
  double sec =  a / globals::tp_1sec; 
  double mins = sec / 60.0;
  double hours = mins / 60.0;
  
  mins -= floor(hours) * 60 ;
  sec  -= floor(hours) * 3600 + floor(mins) * 60;
  
  int h = floor(hours);
  int m = floor(mins);
  int s = floor(sec);

  // return 00:00:00 format
  std::stringstream ss;
  if ( h < 10 ) ss << "0";
  ss << h << ":";
  if ( m < 10 ) ss << "0";
  ss << m << ":";
  if ( s < 10 ) ss << "0";
  ss << s;
  return ss.str();		  
}



std::string Helper::timestring( int h , int m , int s )
{
  // return 00:00:00 format
  std::stringstream ss;
  if ( h < 10 ) ss << "0";
  ss << h << ":";
  if ( m < 10 ) ss << "0";
  ss << m << ":";
  if ( s < 10 ) ss << "0";
  ss << s;
  return ss.str();		  
}


double Helper::position( uint64_t a , uint64_t tot , int * h , int * m , int *s)
{

  // i.e. not that this matters, but for 0-based scaling for 'a'
  --tot;
  
  // a is in tp-units
  double sec = a / (double)globals::tp_1sec;
  double sec2 = sec;
  double mins = sec / 60.0;
  double hours = mins / 60.0;
  
  mins -= floor(hours) * 60 ;
  sec  -= floor(hours) * 3600 + floor(mins) * 60;

  *h = floor(hours);
  *m = floor(mins);
  *s = floor(sec);
  
  return sec2 / (double)(tot/globals::tp_1sec);
}


bool Helper::add_clocktime( int *h , int *m , int *s , uint64_t a , int * msec )
{
  // assumes starting time is always msec == 0
  // (i.e. as per EDF starttime), but this allows for 'a' to imply a fractional time to be output)

  double sec   = *s + (*m)*60 + (*h)*60*60;
  double sec2  = a / (double)globals::tp_1sec;
  double fsec  = sec + sec2;
  
  double fmins  = fsec / 60.0;
  double fhours = fmins / 60.0;
  fmins  -= floor(fhours) * 60 ;
  fsec   -= floor(fhours) * 3600 + floor(fmins) * 60;
  
  // loop around the (24-hr) clock
  if ( fhours > 24 ) fhours -= 24;
  
  *h = floor(fhours);
  *m = floor(fmins);
  *s = floor(fsec);
  
  if ( msec != NULL ) 
    {
      *msec = 1000 * ( fsec - *s );
    }
  
  return true;
}

bool Helper::timestring( const std::string & t, int * h, int *m , int *s)
{
  *h = *m = *s = 0;
  std::vector<std::string> tok = Helper::parse( t , ":.-" );
  if ( tok.size() < 2 ) return false;
  if ( tok.size() > 3 ) return false;
  if ( ! Helper::str2int( tok[0] , h ) ) return false;
  if ( ! Helper::str2int( tok[1] , m ) ) return false;
  if ( tok.size() == 3 )
    if ( ! Helper::str2int( tok[2] , s ) ) return false;  
  return true;
}





void Helper::swap_in_variables( std::string * t , const std::map<std::string,std::string> & vars )
{

  // variable must be in the form   ${var} 
  std::string s;
  for (int i=0;i<t->size();i++)
    {
      if ( (*t)[i] != '$' ) { s = s + (*t)[i]; continue; } 
      ++i;
      if ( i == t->size() ) Helper::halt( "badly formed variable:" + *t );
      if ( (*t)[i] != '{' ) Helper::halt( "badly formed variable:" + *t );
      std::string varname;
      while (1)
	{
	  ++i;
	  if ( i == t->size() ) Helper::halt( "badly formed variable" );
	  
	  if ( (*t)[i] != '}' ) varname += (*t)[i];
	  else
	    {
	      if ( vars.find( varname ) == vars.end() )
		Helper::halt( "variable ${" + varname + "} was not specified" );
	      else // swap in new text
		{
		  s += vars.find( varname )->second;
		  break;
		}
	    }
	}
    }
  *t = s;
}

std::vector<std::string> Helper::file2strvector( const std::string & filename )
{
  if ( ! Helper::fileExists( filename ) ) Helper::halt( "could not find " + filename );
  std::ifstream IN1( filename.c_str() , std::ios::in );
  std::vector<std::string> d;
  while ( ! IN1.eof() )
    {
      std::string x;
      IN1 >> x;
      if ( IN1.eof() ) break;
      d.push_back(x);
    }
  IN1.close();
  return d;
}


bool Helper::hhmmss( const clocktime_t & ct , const interval_t & interval , std::string * t1 , std::string * t2 , const int dp )
{
  *t1 = ".";
  *t2 = ".";
  
  double tp1_sec =  interval.start / (double)globals::tp_1sec;
  clocktime_t present1 = ct;
  present1.advance( tp1_sec / 3600.0 );
  // add down to 1/100th of a second                                                                                                                                               
  double tp1_extra = tp1_sec - (long)tp1_sec;
  
  // rewind stop by 1 unit
  double tp2_sec =  (interval.stop-1LLU) / (double)globals::tp_1sec;
  clocktime_t present2 = ct;
  present2.advance( tp2_sec / 3600.0 );
  double tp2_extra = tp2_sec - (long)tp2_sec;
  
  *t1 = present1.as_string() +  Helper::dbl2str_fixed( tp1_extra , dp  ).substr(1) ;
  *t2 = present2.as_string() +  Helper::dbl2str_fixed( tp2_extra , dp  ).substr(1) ;
  
  return true;
}
