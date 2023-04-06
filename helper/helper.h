
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

#ifndef __HELPER_H__
#define __HELPER_H__

#include <iostream>

#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include <stdint.h>
#include <map>
#include <cmath>

struct interval_t;
struct clocktime_t;
class gzifstream;

namespace Helper 
{

  std::string toupper( const std::string & );  


  // trim from start
  // static inline std::string &ltrim(std::string &s) {
  //   s.erase(s.begin(), std::find_if( s.begin(), s.end(),  std::not1(std::ptr_fun<int, int>(std::isspace))  ));
  //   return s;
  // }

  // trim from end
  // static inline std::string &rtrim( std::string s ) {
  //    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
  //    return s;
  //  }

  // Updated for C++11 using lambda function; trim now do not modify existing string...
  
  // trim from start
  static inline std::string ltrim( std::string s ) {
    s.erase(s.begin(), std::find_if( s.begin(), s.end(),  [](int c) {return !std::isspace(c);} ));
    return s;
  }
 
  // trim from end
  static inline std::string rtrim(std::string s) {
    s.erase(std::find_if( s.rbegin(), s.rend(),  [](int c) {return !std::isspace(c);} ).base(), s.end() );
    return s;
  }
  
  // trim from both ends
  static inline std::string lrtrim( std::string s ) {
    return ltrim(rtrim(s));
  }

  static inline std::string unquote(const std::string &s , const char q2 = '"' ) {
    if ( s.size() == 0 ) return s;
    int a = ( s[0] == '"' || s[0] == q2 ) ? 1 : 0;
    int b = ( s[s.size()-1] == '"' || s[s.size()-1] == q2 ) ? 1 : 0 ;
    return s.substr(a,s.size()-a-b);
  }

  std::string remove_all_quotes(const std::string &s , const char q2 = '"' );
  std::string quote_spaced( const std::string & s );

  std::string quote_if( const std::string & s , char q );
  std::string quote_if( const std::string & s , char q , char p );

  std::string sanitize( const std::string & , const char except );
  std::set<std::string> sanitize( const std::set<std::string> & , const char except  ); 
  std::string sanitize( const std::string & , const std::set<char> * except = NULL );
  std::set<std::string> sanitize( const std::set<std::string> & , const std::set<char> * except = NULL ); 
  
  bool yesno( const std::string & );
  
  std::string search_replace( const std::string & , char a , char b );

  std::string search_replace( std::string s , const std::string & a , const std::string & b , const bool only_root = false );
  
  void swap_in_variables( std::string * , std::map<std::string,std::string> * );

  bool swap_in_includes( std::string * , const std::string & delim = "," );
			
  void expand_numerics( std::string * );

  void process_block_conditionals( std::string * , const std::map<std::string,std::string> & );

  std::string insert_indiv_id( const std::string & id , const std::string & str );

  bool file_extension( const std::string & , const std::string & , bool with_period = true );
  
  bool is_folder( const std::string & f );

  void build_sample_list( const std::vector<std::string> & );

  void merge_EDFs( const std::vector<std::string> & );

  void repath_SL( const std::vector<std::string> & );
  
  void compile_txttabs( const std::string & );

  std::vector<std::string> file2strvector( const std::string & );

  // binary I/O
  
  void bwrite( std::ofstream & O , const std::string & s );
  void bwrite( std::ofstream & O , int i );
  void bwrite( std::ofstream & O , double d );
  std::string bread_str( std::ifstream & I );
  int bread_int( std::ifstream & I );
  double bread_dbl( std::ifstream & I );
  void bskip_dbl( std::ifstream & I , const int n );
  void bskip_int( std::ifstream & I , const int n );

   
  
  // case insenstive string comparison
  bool iequals(const std::string& a, const std::string& b);

  // case insenstive string-root match
  bool imatch(const std::string& a, const std::string& b , unsigned int min = 0 );

  // case-insensitive any match
  bool contains( const std::string& a, const std::string& b );
  
  // print vector
  template <class T> std::string print( const std::vector<T> & x , const std::string & label , const int l )
    {
      std::stringstream ss;
      if ( label != "" )     
	ss << "--- " << label << " ---\n";
      int n = l > x.size() ? x.size() : l ;
      for (int i=0;i<n;i++) ss << i << " [ " << x[i] << " ]\n";
      return ss.str();
    }

  bool fileExists(const std::string &);
  std::string expand( const std::string & f );
  bool deleteFile( const std::string & );
  
  std::vector<std::string> file2strvector( const std::string & filename );

  std::istream& safe_getline(std::istream& is, std::string& t);
  gzifstream & zsafe_getline( gzifstream & is , std::string& t);

  void halt( const std::string & msg );
  void warn( const std::string & msg );
  void debug( const std::string & msg );
  bool realnum(double d);
  bool similar( double a, double b , double EPS = 1e-6 );
  
  std::string int2str(int n);  
  std::string int2str(long n);
  std::string int2str(uint64_t n);
  std::string dbl2str(double n);  
  std::string dbl2str(double n, int dp);  
  std::string dbl2str_fixed(double n, int ch );
  std::string brief( const std::string & , int l = 40);
/*   std::string stringize( const std::vector<std::string> & ); */
/*   std::string stringize( const std::set<std::string> & , const std::string & delim = "," ); */

  std::string readfile( const std::string & file );

  void ascii7( std::string * s , char repl );
  void ascii7( std::vector<char> * s , char repl );
  

  template<typename T> 
    std::string stringize( const T & t , const std::string & delim = "," )
    {
      std::stringstream ss;
      
      typename T::const_iterator tt = t.begin();
      while ( tt != t.end() )
	{
	  if ( tt != t.begin() ) ss << delim;
	  ss << *tt;
	  ++tt;
	}
      return ss.str();
    }
  
  
  std::string trim( const std::string & s , const char c = ' ' , const char d = ' ' );
  std::string format( const std::string & , int indent = 10 , int width = 60, bool no_initial_indent = true );

  std::string squash( const std::string & s , const char c );

  
  bool str2dbl(const std::string & , double * ); 
  bool str2int(const std::string & , int * ); 
  bool str2int64(const std::string & , uint64_t * );
  bool str2signed_int64(const std::string & , int64_t * ); 

  template <class T>
    bool from_string(T& t,
		     const std::string& s,
		     std::ios_base& (*f)(std::ios_base&))
    {
      std::istringstream iss(s);
      return !(iss >> f >> t).fail();
    }
  
  uint64_t sec2tp( double );
  double tp2sec( uint64_t );


  template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
  }

  // vector --> string
  template <class T> std::set<T> vec2set( const std::vector<T> & x )
  {
    std::set<T> s;
    for (int i=0; i<x.size(); i++) s.insert(x[i]);
    return s;
  }

  std::vector<std::string> parse(const std::string & item, const std::string & s = " \t\n" , bool empty = false );
  std::vector<std::string> parse(const std::string & item, const char s , bool empty = false );
  
  std::vector<std::string> quoted_parse(const std::string & item , const std::string & s , const char q = '"' , const char q2 = '\'' , bool empty = false );
  std::vector<std::string> quoted_parse(const std::string & item , const char s , const char q = '"' , const char q2 = '\'' , bool empty = false );

  std::vector<std::string> char_split( const std::string & s , const char c , bool empty );
  std::vector<std::string> char_split( const std::string & s , const char c , const char c2 , bool empty );
  std::vector<std::string> char_split( const std::string & s , const char c , const char c2 , const char c3 , bool empty );

  std::vector<std::string> quoted_char_split( const std::string & s , const char c , const char q , const char q2, bool empty );
  std::vector<std::string> quoted_char_split( const std::string & s , const char c , const char c2 , const char q , const char q2, bool empty );
  std::vector<std::string> quoted_char_split( const std::string & s , const char c , const char c2 , const char c3 , const char q , const char q2, bool empty );


  // time-string
  std::string timestring( uint64_t msec , char delim = '.' , bool fractional = true ); 
  std::string timestring( int h , int m , double s , char delim = '.' , bool fractional = false );
  std::string timestring( const std::string & , const interval_t & , char delim = '.' , const std::string & delim2 = " - " );

  double position( uint64_t a , uint64_t tot , int * h , int * m , double *s);  
  bool timestring( const std::string & , int * h, int *m , double *s );
  bool add_clocktime( int *h , int *m , double *s , uint64_t a );

  // given a clock time and an interval, get an output-friendly timestamp (x2)
  bool hhmmss( const clocktime_t & ct , const interval_t & interval , std::string * t1 , std::string * t2 , const int dp = 4 );


  // pretty print proportion to n dps
  std::string pp( double x );
  
}



struct date_t {

  date_t( const std::string & dt )
  {
    // assumes European dd/mm/yy encoding [ or dd.mm.yy or dd-mm-yy ]
    // also, yy can be yyyy
    // mm can be jan/feb etc or 1/2 etc
    std::vector<std::string> tok = Helper::parse( dt , "./-" );
    if ( tok.size() != 3 ) Helper::halt( "invalid date string: " + dt );

    d=m=y=0;
    
    if ( ! Helper::str2int( tok[0] , &d ) )
      Helper::halt( "invalid day value: " + dt );
    if ( ! Helper::str2int( tok[1] , &m ) )
      {
	std::string mm = Helper::toupper( tok[1] );
	if ( mm.size() == 3 )
	  {
	    if      ( mm == "JAN" ) m = 1;
	    else if ( mm == "FEB" ) m = 2;
	    else if ( mm == "MAR" ) m = 3;
	    else if ( mm == "APR" ) m = 4;
	    else if ( mm == "MAY" ) m = 5;
	    else if ( mm == "JUN" ) m = 6;
	    else if ( mm == "JUL" ) m = 7;
	    else if ( mm == "AUG" ) m = 8;
	    else if ( mm == "SEP" ) m = 9;
	    else if ( mm == "OCT" ) m = 10;
	    else if ( mm == "NOV" ) m = 11;
	    else if ( mm == "DEC" ) m = 12;
	  }
      }

    if ( m == 0 )
      Helper::halt( "invalid month value: " + dt );
    if ( ! Helper::str2int( tok[2] , &y ) )
      Helper::halt( "invalid year value: " + dt );
    
    init();
  }

  int diff( const date_t & rhs ) const;

  // days past 1/1/85
  static int count( const date_t & );
  
  static bool leap_year( const int year )
  {
    return (year % 4 == 0 && year % 100 != 0) || year % 400 == 0 ; 
  }

  static int days_in_month( int mn, int yr )
  {
    static int mlength[] =      { 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
    static int leap_mlength[] = { 0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };    
    return leap_year( yr ) ? leap_mlength[mn] : mlength[mn];    
  }

  // given day count (past 1/1/85) return a printable datestring
  static std::string datestring( int c ); 


  std::string as_string() const
  {
    return Helper::int2str( d ) + "-" + Helper::int2str( m ) + "-" + Helper::int2str( y );
  }
  
  date_t( const int d = 1 , const int m = 1 , const int y = 1985 )
    : d(d) , m(m) , y(y) 
  {
    init();
  }
    
  void init()
  {

    // YY --> YYYY conversion 1985 -- 2084
    //         >=85 --> 19++
    //         < 85 --> 20++
    if      ( y >= 0  && y < 85 ) y += 2000;
    else if ( y >= 85 && y < 100 ) y += 1900;
    
    if ( y < 1985 || y > 3000 )
      Helper::halt( "invalid year (range 1985 - 3000): " + Helper::int2str(y) );
    
    if ( m < 1 || m > 12 )
      Helper::halt( "invalid month (range 1 - 12): " + Helper::int2str(m) );

    if ( d < 1 || d > days_in_month( m , y ) )
      Helper::halt( "invalid day (range 1 - [28-31]): " + Helper::int2str(d) );
  }
   
  int y;
  int m;
  int d;
  
  bool operator<( const date_t & rhs ) const
  {
    if ( y < rhs.y ) return true;
    if ( y > rhs.y ) return false;
    if ( m < rhs.m ) return true;
    if ( m > rhs.m ) return false;
    return d < rhs.d;
  }
  
};



struct clocktime_t
{
  
  // default (midnight)
  // day = days past epoch 1/1/85
  //  missing = 0 (i.e.  == 1/1/85)
  clocktime_t() 
  {
    valid = true;
    d=h=m=0;
    s=0.0;
  }

  bool operator==( const clocktime_t & rhs ) const
  {
    if ( ! valid ) return false;
    if ( ! rhs.valid ) return false;
    if ( d != rhs.d ) return false;
    if ( h != rhs.h ) return false;
    if ( m != rhs.m ) return false;
    if ( fabs( s - rhs.s ) > 1e-12 ) return false;
    return true;
  }
  
  // convert time-string to internal 
  clocktime_t( const std::string & t );

  // convert date-string and time-string to internal
  clocktime_t( const std::string & dt , const std::string & tm );

  void parse_string( const std::string & t ); 
  
  // assume from hours, fractional
  //  clocktime_t( double ); 
  
  clocktime_t( int h, int m, double s ) 
  : valid(true) , d(0) , h(h), m(m), s(s)
  { 
    if ( h < 0 || m < 0 || s < 0 ) valid = false;
    if ( h > 23 || m > 59 || s >= 60.0 ) valid = false;
  } 

  // with day specified
  clocktime_t( int d , int h, int m, double s )
  : valid(true) , d(d) , h(h), m(m), s(s) 
  { 
    if ( d < 0 ) valid = false;
    if ( h < 0 || m < 0 || s < 0 ) valid = false;
    if ( h > 23 || m > 59 || s >= 60.0 ) valid = false;
  } 

  clocktime_t( const clocktime_t & t1 ) { copy(t1); } 

  
  clocktime_t & operator= (const clocktime_t & t1 ) { copy(t1); return *this; }
    
  void copy( const clocktime_t & t1 )
  {
    d = t1.d;
    h = t1.h;
    m = t1.m;
    s = t1.s;
    valid = t1.valid;
  }
  
  void reset() 
  {
    valid = true;
    d = 0; h = 0; m=0; s=0.0;
  }

  bool valid;
  int d;
  int h;
  int m;
  double s;

  std::string as_string( const char tchar = '.' , bool fractional = false ) const
  {
    if ( ! valid ) return "NA";
    return Helper::timestring( h,m,s, tchar , fractional );
  }

  // dd-mm-yyyy-hh:mm:ss
  std::string as_datetime_string( const char tchar = '.' , bool fractional = false ) const
  {
    if ( ! valid ) return "NA";
    return date_t::datestring(d) + "-" + Helper::timestring( h,m,s, tchar , fractional );
  }

  std::string as_numeric_string() const
  {
    if ( ! valid ) return "NA";
    return Helper::dbl2str( hours() );
  }

  // calculate mid-point between two times
  bool midpoint( const clocktime_t & t1 , const clocktime_t & t2 );

  // time as numeric (either in units past epoch (if dr == 0 )
  // or relative to some other reference day ('dr') 
  double minutes( const int dr = 0 ) const;
  double hours( const int dr = 0 ) const;
  double seconds( const int dr = 0 ) const;
  
  int rounded_seconds( const int dr = 0 ) const;
  
  //  void advance_1second();  
  bool convert( double hrs );
  bool convert_seconds( double sec );
  
  void advance_tp( uint64_t tp );
  void advance_hrs( double hrs );
  void advance_seconds( double secs );
  void advance( const clocktime_t & t );
  
  static int earlier( const clocktime_t & t1 , const clocktime_t & t2 );
    
  static double difference_hours( const clocktime_t & t1 , const clocktime_t & t2 );

  static double difference_seconds( const clocktime_t & t1 , const clocktime_t & t2 );

  static double ordered_difference_hours( const clocktime_t & t1 , const clocktime_t & t2 );

  static double ordered_difference_seconds( const clocktime_t & t1 , const clocktime_t & t2 );

  
};




#endif

