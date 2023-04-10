

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
#include "edf/edf.h"
#include "zfstream.h"

#include <cmath>
#include <cstdio>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <streambuf>

#ifndef WINDOWS
#include <wordexp.h>
#endif

extern logger_t logger;

int fn_luna_slbuilder(const char * fpath, const struct stat *ptr, int type );

std::string Helper::toupper( const std::string & s )
{
  std::string j = s;
  for (int i=0;i<j.size();i++) j[i] = std::toupper( s[i] );
  return j;
}

std::string Helper::quote_spaced( const std::string & s )
{
  return quote_if( s , ' ' );
}

std::string Helper::remove_all_quotes(const std::string &s , const char q2 )
{
  const int n = s.size();
  int n2 = 0;
  for (int i=0; i<n; i++) { if ( ! ( s[i] == '"' || s[i] == q2 ) ) ++n2; } 
  if ( n2 == n ) return s;
  std::string r( n2 , ' ' );
  int j = 0;
  for	(int i=0; i<n; i++)
    {
      if ( ! ( s[i] == '"' || s[i] == q2 ) )
	{
	  r[j] = s[i];
	  ++j;
	}
    }
  return r;
}

    
std::string Helper::quote_if( const std::string & s , char q )
{
  // empty strings stay as is
  if ( s == "" ) return s;

  // already quoted?
  if ( s[0] == '"' && s[ s.size() - 1 ] == '"' ) return s;

  // does not contain flagged character: return as is
  if ( s.find( q ) == std::string::npos ) return s;

  // otherwise, place quotes
  return "\"" + s + "\"";
}


std::string Helper::quote_if( const std::string & s , char q , char p )
{
  // empty strings stay as is
  if ( s == "" ) return s;

  // already quoted?
  if ( s[0] == '"' && s[ s.size() - 1 ] == '"' ) return s;

  // does not contain either flagged character: return as is
  if ( s.find( q ) == std::string::npos && s.find( p ) == std::string::npos ) return s;

  // otherwise, place quotes
  return "\"" + s + "\"";
}

std::set<std::string> Helper::sanitize( const std::set<std::string> & s ,
					const char except )
{
  std::set<char> x;
  x.insert( except );
  return Helper::sanitize( s , &x );
}


std::set<std::string> Helper::sanitize( const std::set<std::string> & s ,
					const std::set<char> * except )
{
  std::set<std::string> r;
  std::set<std::string>::const_iterator ss = s.begin();
  while ( ss != s.end() )
    {
      r.insert( Helper::sanitize( *ss , except ) );
      ++ss;
    }
  return r;
}


std::string Helper::sanitize( const std::string & s ,
                              const char except )
{
  std::set<char> x;
  x.insert( except );
  return Helper::sanitize( s , &x );
}


std::string Helper::sanitize( const std::string & s ,
			      const std::set<char> * except )
{

  // nb. does not sanitize comma or pipe or quotes or period
  
  // i.e. and so can be used on expressions that include
  // these operators around labels sig="EEG C3-M2","EEG C4-M1"
  //  --> sig="EEG_C3_M2","EEG_C4_M1"
  
  std::string j = s;
  for (int i=0;i<j.size();i++)
    {
      // an  exception?
      if ( except != NULL && except->find( j[i] ) != except->end() )
	continue;

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

std::string Helper::search_replace( std::string s , const std::string & a , const std::string & b , const bool only_root )
{
  // if only_root == T then only replace if 'a' matches from the start of 's' (i.e. pos == 0 )
 
  
  size_t pos = s.find(a);

  if ( only_root )
    {
      if ( pos == 0 ) s.replace(pos, a.size(), b );
      return s;
    }
  
  while( pos != std::string::npos)
    {
      s.replace(pos, a.size(), b );
      pos = s.find( a, pos + b.size() );
    }

  return s;
}

std::string Helper::expand( const std::string & f )
{

#ifdef WINDOWS
  return f;
#else
  // only expand ~ if first character for home-folder subst
  if ( f.size() == 0 ) return f;
  if ( f[0] != '~' ) return f;
  std::string home = getenv("HOME");
  return home + f.substr(1);
#endif

//   wordexp_t exp_result;
//   // place filepath in quotes to handle spaces, other special chars etc
//   //const int illegal = wordexp( ("\""+f+"\"").c_str(), &exp_result, 0);
//   const int illegal = wordexp( f.c_str(), &exp_result, 0);  
//   if ( illegal )
//     Helper::halt( "code " + Helper::int2str(illegal) + "\nproblem parsing filepath: " + f );  
//   std::string r = exp_result.we_wordv[0];
//   //  for (uint32_t i = 1; i < exp_result.we_wordc; ++i)
//   //    r += " " + std::string(exp_result.we_wordv[i]);
//   wordfree(&exp_result);
//   std::cout << "r [" << r << "]\n";
//   return r;
// #endif
}


bool Helper::is_folder( const std::string & f ) { if ( f.size() == 0 ) return false; return f[f.size()-1]== globals::folder_delimiter; } 


bool Helper::file_extension( const std::string & f, const std::string & ext , bool with_period )
{
  if ( with_period ) 
    {
      int l = ext.size() + 1;  
      if ( f.size() < l ) return false;
      const std::string s = f.substr( f.size() - l );
      return Helper::iequals( s , "." + ext );
    }

  // for matching where extensions and file names differ by a tag,  e.g. file1.edf and file1-annot.xml,  
  // so here match would be "-annot.xml" and would match without a period, rather than "xml" which matches with a period 
  
  else 
    {
      int l = ext.size() ;  
      if ( f.size() < l ) return false;
      const std::string s = f.substr( f.size() - l );
      return Helper::iequals( s , ext );
    }
}




void Helper::halt( const std::string & msg )
{
  
  // some other code handles the exit, e.g. if running under luna-web?
  if ( globals::bail_function != NULL ) 
    globals::bail_function( msg );

  // do not kill the process? i.e. in R mode... likely dangerous...
  if ( ! globals::bail_on_fail ) return;
  
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
  
  return s3.str();
  // std::string retstr = s3.str();
  // return retstr.substr(0,ch);

}

uint64_t Helper::sec2tp( double s )
{

  // to avoid floating point errors, take 's' precision to 1/1000 of a second only, i.e. when reading 
  // input; internally, time-points have 1e-9 precision, so this avoids any floating point issues, 
  // i.e. with small inaccuracies in 's' being scaled up, as double float precision will be better than 1/1000 
  
  if ( s < 0 ) 
    {
      logger << "warning -- cannot have negative time-points, setting to tp=0 (from s=" << Helper::dbl2str( s ) << ")\n";
      return 0; 
    }
  
  int si = floor( s );
  double frac = 1000.0 * ( s - si );
  int ifrac = round( frac * 1000.0 ) / 1000.0;  
  uint64_t retval = si * globals::tp_1sec + ifrac * globals::tp_1000thsec ;   
  //std::cout << std::setprecision(10) << s << "\t" << retval << "\n";
  return retval;
}


double Helper::tp2sec( uint64_t tp )
{
  return (double)tp * globals::tp_duration;
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

bool Helper::str2signed_int64(const std::string & s , int64_t * i)
{
  return from_string<int64_t>(*i,s,std::dec);
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

void Helper::ascii7( std::string * s , char repl )
{
  const int sz = s->size();
  for ( int i=0; i<sz; i++) 
    if ( (*s)[i] < 32 || (*s)[i] > 126 ) (*s)[i] = repl;    
}

void Helper::ascii7( std::vector<char> * s , char repl )
{
  const int sz = s->size();
  for ( int i=0; i<sz; i++) 
    if ( (*s)[i] < 32 || (*s)[i] > 126 ) (*s)[i] = repl;  
}

std::string Helper::squash( const std::string & s , const char c )
{
  std::vector<char> t;
  const int n = s.size();
  for (int i=0; i<n; i++)
    {
      if ( i == 0 ) t.push_back( s[i] );
      else if ( s[i] != c ) t.push_back( s[i] );
      else if ( s[i-1] != c ) t.push_back( s[i] );
    }
  
  std::string ret(t.begin(), t.end());
  return ret;
}

std::string Helper::trim( const std::string & s , const char c , const char d )
{
  int first_nonspace = 0;
  int last_nonspace = s.size()-1;
  for (int i=0;i<s.size();i++) { if ( s[i] == c || s[i] == d ) ++first_nonspace; else break; }
  for (int i=s.size()-1;i!=0;i--) { if ( s[i] == c || s[i] == d ) --last_nonspace; else break; }
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

std::vector<std::string> Helper::parse(const std::string & item, const char s , bool empty )
{
  return Helper::char_split( item , s , empty ); 
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

std::vector<std::string> Helper::quoted_parse(const std::string & item , const char s , const char q , const char q2, bool empty )
{
  return Helper::quoted_char_split( item , s , q, q2, empty ); 
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



// https://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf

std::istream& Helper::safe_getline(std::istream& is, std::string& t)
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


gzifstream & Helper::zsafe_getline( gzifstream & is , std::string& t)
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

      //      std::cout << "c [" << c << "]\n";

      switch (c) 
	{
	case '\n':
	  return is;
	  
	case '\r':
	  if (sb->sgetc() == '\n')
	    sb->sbumpc();
	  return is;


 	case std::streambuf::traits_type::eof() :
	  // replace w/ macro EOF to compile	  
	  // 	case EOF :
 	  // Also handle the case when the last line has no line ending
 	  if(t.empty())
 	    is.setstate(std::ios::eofbit);
 	  return is;
	  
	default:
	  t += (char)c;

	  //	  std::cout << " t [" << t << "]\n";
	}
    }
}


bool Helper::fileExists( const std::string & f )
{

  FILE *file;

  if ( ( file = fopen( f.c_str() , "r" ) ) ) 
    {
      fclose(file);
      return true;
    } 

  return false;

  // std::ifstream inp;
  
  // inp.open(f.c_str(), std::ifstream::in);

  // if(inp.fail())
  //   {
  //     inp.clear(std::ios::failbit);
  //     inp.close();
  //     return false;
  //   }
  // inp.close();
  // return true;

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


bool Helper::imatch(const std::string& a, const std::string& b , unsigned int min )
{
  // only compare up to length of the shortest
  // i.e. 
  //     E
  //     EDF Annotations
  // will match, so  
 
  unsigned int sz = a.size() < b.size() ? a.size() : b.size() ;
  // if specified, require up to 'min' characters
  if ( min != 0 ) sz = min;
  if ( a.size() < min || b.size() < min ) return false;

  for (unsigned int i = 0; i < sz; ++i)
    if (tolower(a[i]) != tolower(b[i]))
      return false;
  return true;
}

bool Helper::yesno( const std::string & s )
{
  // 0 no NO n N F f false FALSE 
  // versus all else  (including empty, i.e. 'var'  --> 'var=T' 
  if ( s.size() == 0 ) return false; // empty == NO
  if ( s[0] == '0' || s[0] == 'n' || s[0] == 'N' || s[0] == 'f' || s[0] == 'F' ) return false;
  return true;
}



std::string date_t::datestring( int c )
{
  // given a count (days past 1/1/85) return a date

  // initiate at c == 0 
  int d = 1;
  int m = 1;
  int y = 1985;

  // add years
  while ( 1 ) {
    int yd = leap_year( y ) ? 366 : 365 ;
    if ( c >= yd )
      {
	c -= yd;
	++y;
      }
    else
      break;
  }

  // add months
  while	( 1 ) {	
    int	md = days_in_month( m , y );
    if ( c >= md )
      {
        c -= md;
	++m;
      }
    else
      break;
  }

  // add days
  d += c;

  return Helper::int2str( d ) + "-" + Helper::int2str( m ) + "-" + Helper::int2str( y );
  
}


int date_t::count( const date_t & dt ) 
{
  
  int days = 0;
  
  // count up until the final year
  for ( int y1 = 1985 ; y1 < dt.y ; y1++ )
    days += leap_year( y1 ) ? 366 : 365 ;

  // count up until the final month
  for ( int m1 = 1 ; m1 < dt.m ; m1++ )
    days += days_in_month( m1 , dt.y );

  // count final days in last month
  days += dt.d;
  
  // is 0-based count, so -1
  // i.e. 1/1/85 == 0 not 1
  return days - 1 ;
    
}
  
int date_t::diff( const date_t & rhs ) const
{
  return count(*this) - count( rhs );
}

clocktime_t::clocktime_t( const std::string & dt ,  const std::string & tm )
{
  // this allow parsing by three chars: . / -
  date_t date( dt );

  const std::string datetime = date.as_string() + "-" + tm;

  // but then we enfore - delimiters, so that the datetime-string reads okay
  parse_string( datetime );
}

clocktime_t::clocktime_t( const std::string & t )
{
  parse_string(t);
}

void clocktime_t::parse_string( const std::string & t )
{
  valid = false;
  // dates? (sep = '/' or '-' only)
  std::vector<std::string> tok = Helper::parse( t , "-/" );
  if ( tok.size() == 1 )
    {
      d=0;
      valid = Helper::timestring( t , &h, &m, &s );
      if ( h < 0 || m < 0 || s < 0 ) valid = false;
      if ( h > 23 || m > 59 || s > 60.0 ) valid = false;
    }
  else if ( tok.size() == 4 )
    {
      date_t dt( tok[0] + "-" + tok[1] + "-" + tok[2] );
      d = date_t::count( dt );
      valid = Helper::timestring( tok[3] , &h, &m, &s );      
      if ( h < 0 || m < 0 || s < 0 ) valid = false;
      if ( h > 23 || m > 59 || s > 60.0 ) valid = false;
    }

}

bool clocktime_t::midpoint( const clocktime_t & t1 , const clocktime_t & t2 )
{
  
  if ( ! ( t1.valid && t2.valid ) )
    {
      valid = false;
      return false;
    }

  // find first
  int e = earlier( t1 , t2 ) ; 
  
  // if t2 actually comes first, use that 
  if ( e == 2 )
    {
      d=t2.d; h=t2.h; m=t2.m; s=t2.s;      
    }
  else // otherwise t1
    {
      d=t1.d; h=t1.h; m=t1.m; s=t1.s;
    }
  
  // time difference (from t1 to t2)
  
  double abs_diff = fabs( difference_hours( t1 , t2 ) ) ;

  // if dayless, then we may need to pick the shorter path
  if ( t1.d == 0 || t2.d == 0 )
    if ( abs_diff > 12 )
      abs_diff = 24.0 - abs_diff ; 
  
  // advance from the earlier point
  advance_hrs( abs_diff / 2.0 );
  
  return true;
}


double clocktime_t::minutes( const int dr ) const 
{
  return (d-dr)*24*60 + h*60 + m + s/(double)60; 
}

double clocktime_t::hours( const int dr ) const
{
  return (d-dr)*24 + h + m/(double)60 + s/(double)(60*60);
}

double clocktime_t::seconds( const int dr ) const
{
  return (d-dr)*24*60*60 + h*60*60 + m*60 + s ;
}

int clocktime_t::rounded_seconds( const int dr ) const
{
  int si = floor(s);
  if ( s - si > 0.5 ) ++si;
  return (d-dr)*24*60*60 + h*60*60 + m*60 + si ;
}

// convert only time of current day, i.e. ignores day
bool clocktime_t::convert( double hrs )
{
  valid = true;
  if ( hrs < 0 ) valid = false;
  if ( hrs > 24 ) valid = false;
  if ( ! valid ) return false;
  
  double t_hours = hrs;
  double t_mins  = hrs * 60.0;
  double t_secs  = hrs * 3600.0;
    
  t_mins -= floor(t_hours) * 60 ;
  t_secs -= floor(t_hours) * 3600 + floor(t_mins) * 60;
  
  h = floor(t_hours);
  m = floor(t_mins);
  s = t_secs;
  
  return true;
  
}
  
// convert only time of current day (i.e. ignores day)
bool clocktime_t::convert_seconds( double sec ) 
{
  valid = true;
  if ( sec < 0 ) valid = false;
  if ( sec > 86400 ) valid = false;
  if ( ! valid ) return false;
  
  double t_hours = floor( sec / 3600.0 );
  sec -= t_hours * 3600.0;
  
  double t_mins  = floor( sec / 60.0 );
  sec -= t_mins * 60.0;
  
  h = t_hours;
  m = t_mins;
  s = sec;
  
  return true;
  
}


void clocktime_t::advance_tp( uint64_t tp )
{
  // convert to hours
  double sec = tp / globals::tp_1sec;
  advance_seconds( sec );
}


void clocktime_t::advance_hrs( double hrs ) 
{
  if ( hrs == 0 ) return;
      
  // anchor on current day
  double t_hrs = hours( d ); 
  
  t_hrs += hrs; 
  
  // need to wrap around days?
  // adjusting 'd' if appropriate 
  while ( 1 ) 
    { 
      if ( t_hrs >= 0 && t_hrs < 24 ) break; 
      if ( t_hrs < 0 ) { t_hrs += 24.0; if ( d ) --d; }
      else if ( t_hrs >= 24 ) { t_hrs -= 24.0; if ( d ) ++d; } 
    }     
  
  // update this time back to usual format 
  convert( t_hrs ); 
} 

void clocktime_t::advance_seconds( double secs )
{
  
  // use current day as the anchor
  double t_sec = seconds( d );
  
  t_sec += secs;
  
  // need to wrap? (86400 seconds in a day)
  // (adjusting 'd' if appropriate)
  while ( 1 )
    {
      
      if ( t_sec >= 0 && t_sec < 86400 )
	break;
      
      if ( t_sec < 0 )
	{
	  t_sec += 86400.0;
	  if ( d ) --d;
	} 
      else if ( t_sec >= 86400.0 )
	{
	  t_sec -= 86400.0;
	  if ( d ) ++d;
	} 
    }    
  
  // update this time back to usual format
  convert_seconds( t_sec );
}


void clocktime_t::advance( const clocktime_t & t ) 
{
  
  if ( ! t.valid ) { valid = false; return; }

  // the second clock time must note contain a date (i.e. this is an interval
  // of time to advance, not a date per se);  also, it should not be negative
  
  if ( t.d )
    Helper::halt( "internal error: cannot specify a clocktime with a date as the durtion of an event, i.e. need h:m:s format)" );
  
  // get only seconds past the start of this day (i.e. by adding 'd' option)  
  double secs = seconds( d );
  double secs2 = t.seconds(); // this should not have any day specified (see above)

  if ( secs2 < 0 )
    Helper::halt( "internal error: clocktime_t::advance() expects positive values only" );
  
  // advance
  secs += secs2;
  
  // need to wrap? 
  // seconds in the day = 0 .. 86400
  // and move day forward (if we are tracking d)
  
  while ( 1 )
    {
      if ( secs >= 86400 )
	{
	  secs -= 86400;
	  if (d) ++d;
	} 
      else break;
    }
  
  // convert back to h/m/s
  
  double t_secs   = secs;
  double t_mins   = secs / 60.0;
  double t_hours  = secs / 3600.0;
  
  t_mins -= floor(t_hours) * 60 ;
  t_secs -= floor(t_hours) * 3600 + floor(t_mins) * 60;
  
  h = floor(t_hours);
  m = floor(t_mins);
  s = t_secs;
  
}

int clocktime_t::earlier( const clocktime_t & t1 , const clocktime_t & t2 )
{

  // exact match?
  if ( t1 == t2 )
    {
      return 0;
    }

  // if day is specified for both times, then this comparison is unambiguous

  // otherwise, if day is not known, the we order the two times to give the shortest
  // difference, i.e.:

  //   t1 = 09:00  t2 = 09:10    [ t1 is earlier,  as 10 mins < 23hours+50mins ]
  //   t1 = 22:00  t2 = 02:00    [ t1 is earlier as 4 hrs < 20 hrs ] 
  
  const bool dayless = t1.d == 0 || t2.d == 0 ;
  
  // unambiguous based on day differences?
  if ( ! dayless )
    {
      //      std::cout << " day comp " << t1.d << " " << t2.d << "\n";
      if ( t1.d < t2.d ) { return 1; }
      if ( t2.d < t1.d ) { return 2; } 
      //std::cout << " still not dayless\n";
      
      // if one the same day, and day is explicitly specified, then difference_hour() sign
      // is also unambiguous

      double d1 = difference_hours( t1 , t2 ) ;
      
      return d1 < 0 ? 2 : 1 ;
      
    }

  // otherwise, if we are dayless, then we need to resolve e.g. that 23:00 comes "before' 01:00
  // using the shortest-duration rule;
  
  // i.e. flip is |diff| is > 12 based on the t2 - t1 difference
  
  double d1 = difference_hours( t1 , t2 )  ;
  double dabs = fabs( d1 );
  //std::cout << "d1 dabs = " << d1 << " " << dabs << "\n";
  if ( dabs <= 12 )
    return d1 < 0 ? 2 : 1 ;
  else
    return d1 < 0 ? 1 : 2 ;

  return 0;
}


double clocktime_t::difference_hours( const clocktime_t & t1 , const clocktime_t & t2 )
{

  // if days specified, use t1 as anchor

  const bool dayless = t1.d == 0 || t2.d == 0 ;

  // if no day informtion, then
  double t1h = 0 , t2h = 0;

  if ( dayless )
    {
      // in case one of the days has a day code, we should supply
      // both explicitly
      t1h = t1.hours( t1.d );
      t2h = t2.hours( t2.d );
    }
  else
    {
      // when both specified, t1 will be the anchor day
      t1h = t1.hours( t1.d );
      t2h = t2.hours( t1.d );
    }
  
  //  std::cout << "difference_hours(): t1 = " << t1h << " t2 = " << t2h << " --> " << t2h - t1h << "\n";
  
  // signed difference
  return t2h - t1h;
  
}


double clocktime_t::ordered_difference_hours( const clocktime_t & t1 , const clocktime_t & t2 )
{

  // exact match?
  if ( t1 == t2 ) return 0;

  // here we *assume* t1 happens before t2
  
  // thus  22 8   means from 22 to 8 (not 8 to 22)
  // and assume time wrap is always within 1 day

  double t1h = t1.hours();
  double t2h = t2.hours();
  
  if ( t2h < t1h ) // e.g. 22 to 8 , means it wraps
    {
      return 24.0 - t1h + t2h;
    }
  else // no overnight wrap, 
    {
      return t2h - t1h;
    }
  
  return 0;
  
}

double clocktime_t::ordered_difference_seconds( const clocktime_t & t1 , const clocktime_t & t2 )
{

  if ( t1 == t2 )
    {
      //std::cout << " IS EQ\n";
      return 0;
    }
  
  // here we *assume* t1 happens before t2
  // thus  22 8   means from 22 to 8 (not 8 to 22)
  // and assume time wrap is always within 1 day
  // i.e. no use of date info for now
  
  double t1s = t1.seconds();
  double t2s = t2.seconds();

  //  std::cout << " t1, t2 = " << t1s << " " << t2s << "\n";
  
  if ( t2s < t1s ) 
    {
      return 86400.0 - t1s + t2s;
    }
  else // no overnight wrap, 
    {
      return t2s - t1s;
    }
  
  return 0;
  
}


double clocktime_t::difference_seconds( const clocktime_t & t1 , const clocktime_t & t2 )
{

  // if days specified, use t1 as anchor

  const bool dayless = t1.d == 0 || t2.d == 0 ;

  // if no day informtion, then
  double t1s = 0 , t2s = 0;

  if ( dayless )
    {
      // in case one of the days has a day code, we should supply
      // both explicitly
      t1s = t1.seconds( t1.d );
      t2s = t2.seconds( t2.d );
    }
  else
    {
      // when both specified, t1 will be the anchor day
      t1s = t1.seconds( t1.d );
      t2s = t2.seconds( t1.d );
    }
  
  // signed difference
  return t2s - t1s;

}



std::string Helper::timestring( const std::string & st , 
				const interval_t & i , 
				char delim , 
				const std::string & delim2 )
{  
  int h0 = 0, m0 = 0;
  double s0 = 0.0;
  if ( ! Helper::timestring( st, &h0, &m0, &s0 ) ) return ".";  
  int h1 = h0, m1 = m0; double s1 = s0;
  int h2 = h0, m2 = m0; double s2 = s0;
  Helper::add_clocktime( &h1, &m1, &s1 , i.start );
  Helper::add_clocktime( &h2, &m2, &s2 , i.stop ); 
  
  std::stringstream str;
  str.precision(0);
  str << std::fixed;
  str << Helper::timestring(h1,m1,s1,delim) << delim2 << Helper::timestring(h2,m2,s2,delim);
  return str.str(); 
}


std::string Helper::timestring( uint64_t a , char delim , bool fractional )
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

  // return 00:00:00 or 00:00:00.00 format
  std::stringstream ss;
  if ( h < 10 ) ss << "0";
  ss << h << delim;

  if ( m < 10 ) ss << "0";
  ss << m << delim;

  if ( s < 10.0 ) ss << "0";
  if ( fractional ) 
    ss << std::fixed << std::setprecision( globals::time_format_dp ) << sec;
  else
    ss << s;

  return ss.str();		  
}



std::string Helper::timestring( int h , int m , double sec , char delim , bool fractional )
{
  // adjust any small rounding error
  if ( sec < 0 ) sec = 0;
  
  // return 00:00:00 or 00:00:00.00 format
  std::stringstream ss;

  if ( h < 10 ) ss << "0";
  ss << h << delim;
  if ( m < 10 ) ss << "0";
  ss << m << delim;
  if ( sec < 10.0 ) ss << "0";
  
  if ( fractional ) 
    ss << std::fixed << std::setprecision( globals::time_format_dp ) << sec;
  else
    ss << floor( sec );
  
  return ss.str();		  
}


double Helper::position( uint64_t a , uint64_t tot , int * h , int * m , double *s)
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
  *s = sec;
  
  return sec2 / (double)(tot/globals::tp_1sec);
}


bool Helper::add_clocktime( int *h , int *m , double *s , uint64_t a  )
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
  *s = fsec;
    
  return true;
}

bool Helper::timestring( const std::string & t0, int * h, int *m , double *s )
{
  // primary function to parse a time-string (text) to internal format

  *h = *m = 0;  
  *s = 0.0;

  // because of EDF spec, '.' a valid delimiter for hh.mm.ss
  // makes fractional seconds hard to parse, if we also allow days
  // so, parse separately: if 1+ colon, assume consistent colon-formatting 
  
  // valid formats:     hh:mm          hh.mm
  //                    hh:mm:ss       hh.mm.ss
  //                    hh:mm:ss.ss    hh.mm.ss.ss
  //                 dd:hh:mm:ss.ss   

  //
  // allow AM/PM modifiers, otherwise default assumption is the 24-hour clock
  // (both upper and lower case)
  //
  
  bool am = t0.find( "AM" ) != std::string::npos || t0.find( "am" ) != std::string::npos ;
  bool pm = t0.find( "PM" ) != std::string::npos || t0.find( "pm" ) != std::string::npos ;
  if ( am && pm ) return false;

  std::string t = ( am || pm ) ? "" : t0;

  // strip AM/PM and any spaces
  if ( am || pm )
    {
      for (int c=0;c<t0.size();c++)
	{
	  if ( t0[c] == ' ' ) continue;
	  if ( t0[c] == 'P' || t0[c] == 'A' || t0[c] == 'M' ) continue;
	  if ( t0[c] == 'p' || t0[c] == 'a' || t0[c] == 'm' ) continue;
	  t += t0[c];	  
	}
    }
  
  //
  // is this colon-delimited, or period-delimited?
  //
  
  std::vector<std::string> tokc = Helper::parse( t , ":" );

  // colon-delimited?
  if ( tokc.size() > 1 ) 
    {
      if ( tokc.size() == 2 ) // hh:mm
	{	  
	  if ( ! Helper::str2int( tokc[0] , h ) ) return false;
	  if ( ! Helper::str2int( tokc[1] , m ) ) return false;
	  if ( am || pm )
	    {
	      if ( *h < 1 || *h > 12 ) return false;
	      if ( pm ) *h += 12;
	      if ( *h == 24 ) *h = 0;
	    }
	  return true;
	}
      else if ( tokc.size() == 3 ) // hh:mm:ss
	{
	  if ( ! Helper::str2int( tokc[0] , h ) ) return false;
	  if ( ! Helper::str2int( tokc[1] , m ) ) return false;
	  if ( ! Helper::str2dbl( tokc[2] , s ) ) return false;
	  if ( am || pm )
	    {
	      if ( *h < 1 || *h > 12 ) return false;
	      if ( pm ) *h += 12;
	      if ( *h == 24 ) *h = 0;
	    }
	  return true;
	}
      else if ( tokc.size() == 4 ) // dd:hh:mm:ss
	{
	  int day = 0;
	  if ( ! Helper::str2int( tokc[0] , &day ) ) return false;
	  if ( ! Helper::str2int( tokc[1] , h ) ) return false;
	  if ( ! Helper::str2int( tokc[2] , m ) ) return false;
	  if ( ! Helper::str2dbl( tokc[3] , s ) ) return false;

	  if ( am || pm )
	    {
	      if ( *h < 1 || *h > 12 ) return false;
	      if ( pm ) *h += 12;
	      if ( *h == 24 ) *h = 0;
	    }

	  // adjust by day
	  h += 24 * day;
	  return true;
	}
      else 
	return false;      
    }
  

  //
  // if no colons, assume this is '.' (or '-') delimited
  //
  
  std::vector<std::string> tok = Helper::parse( t , ".-" );

  if ( tok.size() == 2  ) // hh.mm
    {
      if ( ! Helper::str2int( tok[0] , h ) ) return false;
      if ( ! Helper::str2int( tok[1] , m ) ) return false;
      if ( am || pm )
	{
	  if ( *h < 1 || *h > 12 ) return false;
	  if ( pm ) *h += 12;
	  if ( *h == 24 ) *h = 0;
	}
      return true;
    }
  else if ( tok.size() == 3 ) // hh.mm.ss
    {      
      if ( ! Helper::str2int( tok[0] , h ) ) return false;
      if ( ! Helper::str2int( tok[1] , m ) ) return false;
      if ( ! Helper::str2dbl( tok[2] , s ) ) return false;        
      if ( am || pm )
	{
	  if ( *h < 1 || *h > 12 ) return false;
	  if ( pm ) *h += 12;
	  if ( *h == 24 ) *h = 0;
	}      
      return true;
    }
  else if ( tok.size() == 4 ) // hh.mm.ss.ss
    {
      if ( ! Helper::str2int( tok[0] , h ) ) return false;
      if ( ! Helper::str2int( tok[1] , m ) ) return false;
      if ( ! Helper::str2dbl( tok[2] + "." + tok[3] , s ) ) return false;
      if ( am || pm )
	{
	  if ( *h < 1 || *h > 12 ) return false;
	  if ( pm ) *h += 12;
	  if ( *h == 24 ) *h = 0;
	}
      return true;
    }
  
  return false;

}


std::string Helper::insert_indiv_id( const std::string & id , const std::string & str )
{
  // replace all instances of 'globals::indiv_wildcard' in 'str' with 'id'  
  std::string v = str;
  while ( v.find( globals::indiv_wildcard ) != std::string::npos )
    {
      int p = v.find( globals::indiv_wildcard  );
      v = v.substr( 0 , p ) + id + v.substr(p+1);
    }
  return v;
}


void Helper::expand_numerics( std::string * t )
{
  // expand [SIG][1:4] to SIG1,SIG2,SIG3,SIG4  
  // search for '][' 
  
  std::map<int,int> splices; // txt locations
  // maps keyed on start position of 'splices'
  std::map<int,int> starts, stops; // to-be-expanded numbers, e.g. 1 and 4
  std::map<int,std::string> root; // e.g. SIG
  
  for (int i=1;i<t->size();i++)
    {
      if ( (*t)[i-1] == ']' && (*t)[i] == '[' ) 
	{
	  int j=i-1;
	  while ( 1 ) { 
	    --j;
	    if ( j < 0 ) Helper::halt( "bad format for [var][n:m]" );
	    if ( (*t)[j] == '[' ) break;
	  }
	  
	  int k=i+1;
	  while ( 1 ) { 
	    ++k;
	    if ( k == t->size() ) Helper::halt( "bad format for [var][n:m]" );
	    if ( (*t)[k] == ']' ) break;
	  }
	  
	  std::string s = t->substr( j , k - j + 1 );
	  std::vector<std::string> tok = Helper::parse( s , "][" );
	  if ( tok.size() != 2 ) Helper::halt( "bad format for [var][n:m]" );
	  
	  std::vector<std::string> tok2 = Helper::parse( tok[1] , ":" );
	  if ( tok2.size() != 2 ) Helper::halt( "bad format for [var][n:m]" );
	  int s1 , s2;
	  if ( ! Helper::str2int( tok2[0] , &s1 ) ) Helper::halt( "bad format for [var][n:m]" );
	  if ( ! Helper::str2int( tok2[1] , &s2 ) ) Helper::halt( "bad format for [var][n:m]" );
	  starts[j] = s1;
	  stops[j] = s2;
	  root[j] = tok[0];
	  splices[j] = k;
	}
    }

//   std::map<int,int>::const_iterator jj = starts.begin();
//   while ( jj != starts.end() )
//     {
//       std::cout << "--> " 
// 		<< jj->second << "\t"
// 		<< stops[ jj->first ] << "\t"
// 		<< root[ jj->first ] << "\n";
//       ++jj;
//     }

  //
  // splice in...
  //
  
  if ( root.size() == 0 ) return;
  
  int p = 0;
  std::string s;
  std::map<int,int>::const_iterator jj = splices.begin();
  while ( jj != splices.end() )
    {
      s += t->substr( p , jj->first - p ); 
      int d = starts[ jj->first ] < stops[ jj->first ] ? 1 : -1 ; 
      bool first = true;
      for (int a = starts[ jj->first ] ; a <= stops[ jj->first ] ; a += d ) 
	{
	  if ( ! first ) s += ",";
	  s += root[ jj->first ] + Helper::int2str( a );
	  first = false;
	}
      // skip to end 
      p = splices[ jj->first ] + 1;
      ++jj;
    }
  // now add final part
  s += t->substr( p );
  
  // std::cout << "s = [" << s << "]\n";
  // std::cout << "t = [" << *t << "]\n";
  *t = s ; 
}

void Helper::swap_in_variables( std::string * t , std::map<std::string,std::string> * vars )
{

  // variable must be in the form   ${var} 
  // definitions can be as ${var=values,etc}
    
  int open = 0;
  std::string s;
  for (int i=0;i<t->size();i++)
    {
      if ( (*t)[i] != '$' ) { s = s + (*t)[i]; continue; } 
      ++i;
      
      if ( i == t->size() ) Helper::halt( "badly formed variable:" + *t );
      if ( (*t)[i] != '{' ) Helper::halt( "badly formed variable:" + *t );
      ++open;
      std::string varname;
      while (1)
	{
	  ++i;
	  if ( i == t->size() ) Helper::halt( "badly formed variable" );
	  	  
	  if ( (*t)[i] != '}' || open > 1 ) 
	    {
	      varname += (*t)[i];
	      if ( (*t)[i] == '}' ) --open;
	      if ( (*t)[i] == '{' ) ++open;
	    }	  
	  else
	    {	      
	      open = 0;

	      // special?
	      if ( cmd_t::is_special( varname ) )
		{
		  Helper::halt( varname + " is a reserved variable and cannot be used in a script" );
		}
	      // definition?
	      else if ( varname.find( "=" ) != std::string::npos )
		{
		  std::vector<std::string> tok = Helper::parse( varname , "=" );
		  
		  if ( tok.size() != 2 ) Helper::halt( "bad format for ${var=value} definition" );
		  // recursively swap in any existing variables in the defiinition
		  // ${a=${b}} 
		  Helper::swap_in_variables( &tok[1] , vars );
		  (*vars)[ tok[0] ] = tok[1]; 
		  break;
		}
  	      else if ( vars->find( varname ) == vars->end() )
		Helper::halt( "variable ${" + varname + "} was not specified" );
	      else // swap in new text
		{
		  s += vars->find( varname )->second;
		  break;
		}
	    }
	}
    }
  *t = s;

}


bool Helper::swap_in_includes( std::string * t ,			       
			       const std::string & delim )
{

  bool changed = false;
  
  // includes must be in the form @{include} 

  std::string s;

  for (int i=0;i<t->size();i++)
    {
      
      if ( (*t)[i] != '@' ) { s = s + (*t)[i]; continue; }       
      ++i;
      changed = true;

      if ( i == t->size() ) Helper::halt( "badly formed @{include}:" + *t );
      if ( (*t)[i] != '{' ) Helper::halt( "badly formed @{include}:" + *t );
      
      std::string filename;
      while (1)
	{
	  ++i;
	  if ( i == t->size() ) Helper::halt( "badly formed @{include}" );
	  	  
	  if ( (*t)[i] != '}' ) filename += (*t)[i];
	  else break;	  
	}
      
      // check for inserting file contents
      if ( ! Helper::fileExists( filename ) )
	Helper::halt( "could not find @{include} file: " + filename );
      
      std::string insert;
      std::ifstream IN( filename.c_str() , std::ios::in );
      while ( ! IN.eof() )
	{
	  std::string item;
	  IN >> item;
	  if ( IN.eof() ) break;
	  if ( insert != "" ) insert += delim ;
	  insert += item;
	}
      IN.close();
      s += insert;

      // continue on to the next character
    }

  // all done
  *t = s;

  return changed;
}
  

void Helper::process_block_conditionals( std::string * t , const std::map<std::string,std::string> & vars )
{
  
  // either tag=1 or tag=0 and tag2=1 etc
  // or  add=tag,tag2
  
  std::set<std::string> adds;
  if ( vars.find( "add" ) != vars.end() )
    {
      std::vector<std::string> tok = Helper::parse( vars.find( "add" )->second , "," );
      for (int i=0; i<tok.size(); i++) {
	if ( cmd_t::is_special( tok[i] ) ) Helper::halt( "cannot specify special variable " + tok[i] );
	adds.insert( tok[i] );
      }

    }

  // 
  // [[var 
  //    ...
  //    include ...
  //    ...
  // ]]var
  //
    
  
  std::string s;
  bool include = true;

  std::set<std::string> includes;
  std::set<std::string> excludes;

  for (int i=0;i<t->size();i++)
    {
      
      //
      // end of an inclusion block?
      //
      
      if ( i < t->size() - 1 && (*t)[i] == ']' && (*t)[i+1] == ']' ) 
	{
	  ++i;
	  std::string h = "";
	  while (1) 
	    { 
	      ++i;
	      if ( i == t->size() ) break; // can be EOF
	      if ( (*t)[i] == ' ' || (*t)[i] == '\t' || (*t)[i] == '\n' ) break;
	      h += (*t)[i];
	    }	  

	  if ( cmd_t::is_special( h ) ) 
	    Helper::halt( h + " is a special reserved variable, cannot be used for a block-conditional" );
	  
	  bool was_exclude = excludes.find(h) != excludes.end();
	  bool was_include = includes.find(h) != includes.end();
	  
	  if ( was_exclude ) excludes.erase( excludes.find(h)  );
	  else if ( was_include ) includes.erase( includes.find(h) );
	  
	  // set current inclusion status
	  include = excludes.size() == 0 ;

	  continue;
	}

      
      //
      // skipping?
      //

      if ( ! include ) continue;

      //
      // start of inclusion block?
      //

      if ( include && i < t->size() - 1 && (*t)[i] == '[' && (*t)[i+1] == '[' ) 
	{
	  ++i;
	  std::string h = "";
	  while (1) 
	    { 
	      ++i;
	      if ( i == t->size() ) Helper::halt( "badly formed inclusion block" );
	      if ( (*t)[i] == ' ' || (*t)[i] == '\t' || (*t)[i] == '\n' ) break;
	      h += (*t)[i];
	    }
	  
	  bool v_inc = vars.find( h ) != vars.end() && vars.find( h )->second != "0";
	  bool a_inc = adds.find( h ) != adds.end();

	  include = v_inc || a_inc;

	  if ( includes.find( h ) != includes.end() || excludes.find(h) != excludes.end() )
	    Helper::halt( "bad format for conditional block: [["+h + " already set" );
	  
	  if ( include ) 
	    includes.insert( h );
	  else 
	    excludes.insert( h );


	  continue;
	}
      
      
      //
      // normal add
      //

      s += (*t)[i];
      
    }

  // copy back
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


bool Helper::hhmmss( const clocktime_t & ct , 
		     const interval_t & interval , 
		     std::string * t1 , std::string * t2 , 
		     const int dp )
{
  *t1 = ".";
  *t2 = ".";
  
  double tp1_sec =  interval.start / (double)globals::tp_1sec;
  clocktime_t present1 = ct;
  present1.advance_seconds( tp1_sec );
  // add down to 1/100th of a second                                                                                                                                               
  double tp1_extra = tp1_sec - (long)tp1_sec;
  
  // rewind stop by 1 unit
  double tp2_sec =  (interval.stop-1LLU) / (double)globals::tp_1sec;
  clocktime_t present2 = ct;
  present2.advance_seconds( tp2_sec );
  double tp2_extra = tp2_sec - (long)tp2_sec;
  
  *t1 = present1.as_string(':') +  Helper::dbl2str_fixed( tp1_extra , dp  ).substr(1) ;
  *t2 = present2.as_string(':') +  Helper::dbl2str_fixed( tp2_extra , dp  ).substr(1) ;
  
  return true;
}




void Helper::compile_txttabs( const std::string & d )
{
  // compile all files across all subfolders of d 
  // make a special folder '00_all'
  
  // create folder if it does not exist
  std::string syscmd = globals::mkdir_command + " " + d + "/00_all";   
  int retval = system( syscmd.c_str() );
  
  // 1) get all subfolders (indivs)
  
  

}


bool Helper::contains( const std::string & a , const std::string & b )
{
  std::string au = Helper::toupper( a );
  std::string bu = Helper::toupper( b );
  return au.find( bu ) != std::string::npos;
}





void Helper::repath_SL( const std::vector<std::string> & tok )
{
  // A -- B and use std::cin for SL
  if ( tok.size() != 2 )
    Helper::halt( "expecting exactly two arguments: old-path new-path < s.lst > new.lst" );
  
  const std::string s1 = tok[0];
  const std::string s2 = tok[1];

  while ( 1 )
    {
      std::string line;
      Helper::safe_getline( std::cin , line );
      if ( std::cin.bad() || std::cin.eof() ) break;
      if ( line == "" ) continue;

      // expect 2+ cols
      std::vector<std::string> tok2 = Helper::parse( line , "\t" );
      if ( tok2.size() < 2 )
	Helper::halt( "requires (ID) | EDF file | (optional ANNOT files)" );
      
      // if s1 is '.' then means just append this prefix (unless starts /)
      if ( s1 == "." )
	{
	  // ensure we have a folder delim at end of s2
	  std::string xdelim = "";
	  if ( s2[ s2.size() - 1 ] != globals::folder_delimiter )
	    xdelim = globals::folder_delimiter ;
	  
	  for (int i=1; i<tok2.size(); i++)
            if ( tok2[i][0] != globals::folder_delimiter )
	      tok2[i] = s2 + xdelim + tok2[i];
	}
      else // replace cols 2+ ( true in search_replace() implies we only match roots of strings)
	{
	  for (int i=1; i<tok2.size(); i++)
	    tok2[i] = Helper::search_replace( tok2[i] , s1 , s2 , true );
	}
      
      for (int i=0; i<tok2.size(); i++)
	std::cout << ( i != 0 ? "\t" : "" ) << tok2[i] ;
      std::cout << "\n";
    }

  return;
}

std::string Helper::readfile( const std::string & f )
{
  std::string s;

  std::string filename = Helper::expand( f );

  if ( ! Helper::fileExists( filename ) ) return s;

  std::ifstream IN1( f.c_str() , std::ios::in );

  while ( 1 )
    {
      std::string line;
      Helper::safe_getline( IN1 , line );
      if ( IN1.bad() || IN1.eof() ) break;
      if ( line == "" ) continue;

      s += line + "\n";
    }
  IN1.close();
  
  return s;
  
}


std::string Helper::pp( double x )
{
  // print a proportion 0 .. 1  w/ only 2 dps
  int x2 = round( 100 * x );
  return Helper::dbl2str( x2 / 100.0 );
}





void Helper::bwrite( std::ofstream & O , const std::string & s ) 
{
  uint8_t l = s.size();
  O.write( (char*)( &l ), sizeof(uint8_t) );
  O.write( s.c_str(), l );
}

void Helper::bwrite( std::ofstream & O , int i ) 
{
  O.write( (char*)( &i ), sizeof(int) );
}

void Helper::bwrite( std::ofstream & O , double d ) 
{
  O.write( (char*)( &d ), sizeof(double) );
}

std::string Helper::bread_str( std::ifstream & I )
{
  uint8_t len;
  I.read( (char*)( &len ), sizeof(uint8_t) );
  std::vector<char> b( len );
  I.read( &b[0] , len );
  std::string s( b.begin() , b.end() );
  return s;
}

int Helper::bread_int( std::ifstream & I )
{
  int i;
  I.read( (char*)( &i ), sizeof(int) );
  return i;
}

double Helper::bread_dbl( std::ifstream & I )
{
  double d;
  I.read( (char*)( &d ), sizeof(double) );
  return d;
}

void Helper::bskip_dbl( std::ifstream & I , const int n )
{
  std::vector<double> dummy( n ) ;
  I.read( (char*)( &dummy[0] ), n * sizeof(double) );
}

void Helper::bskip_int( std::ifstream & I , const int n )
{
  std::vector<double> dummy( n ) ;
  I.read( (char*)( &dummy[0] ), n * sizeof(int) );
}



