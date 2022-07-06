
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
#include <vector>
#include <cstring>
#include <fstream>


std::istream& safe_getline(std::istream& is, std::string& t);

std::vector<std::string> char_split( const std::string & s , const char c , bool empty = true );

struct reg_t {  

  reg_t( double a, double b) : start(a), stop(b) { }
  
  double start;
  double stop;

  bool operator<( const reg_t & rhs ) const 
  {
    if ( start < rhs.start ) return true;
    if ( start > rhs.start ) return false;
    return stop < rhs.stop;
  }
};

bool region_is_included( const reg_t & r , const std::set<reg_t> & d );
std::set<reg_t> flatten( const std::set<reg_t> & x );
  


int main(int argc , char ** argv )
{
  
  // regional list-of-regions col1 col2 COL < input.txt > output.txt 

  // if COL == -X means exclude row
  //           -I mean include row
  // else add COL (in header) called this, w/ 0/1 values 
  
  // output:
  //  input with extra column added  (or rows deleted)
  
  if ( argc != 5 )
    {
      std::cerr << " error: usage \n"
		<< "   regional <regions> <start-col> <stop-col> <new-col> < input > output\n";
      std::exit(1);
    }

  // expect just start / stop 
  const std::string region_file = argv[1] ;

  // start (in input data) [ base 1 ] 
  const int col1 = atoi( argv[2] );
  
  // stop (in input data) [ base 1 ] 
  const int col2 = atoi( argv[3] );
  
  const std::string newcol = argv[4];

  const bool exclude_matches = newcol == "-X" ;

  const	bool include_matches = newcol == "-I" ;

  const bool append_col = ! ( exclude_matches || include_matches );
  
  //
  // read regions
  //

  std::set<reg_t> regions;

  std::ifstream IN( region_file.c_str() , std::ios::in );
  while ( 1 )
    {
      double x,y;
      IN >> x >> y;
      if ( IN.eof() || IN.bad() ) break;
      regions.insert( reg_t( x, y ) );
    }

  std::cerr << "read " << regions.size() << " regions\n";

  regions = flatten( regions );

  std::cerr << "after flattening, " << regions.size() << " regions\n";
    
  //
  // Read data 
  //

  int data_line = 0 , included_data_line = 0;
  
  while ( 1 )
    {

      // read line
      std::string line;
      safe_getline( std::cin , line );
      if ( std::cin.eof() ) break;
      std::vector<std::string> tok = char_split( line , '\t' );
      if ( line == "" || tok.size() == 0 ) continue;
      
      // handle header/comment lines in common files
      const bool header = tok[0] == "ID" ;
      const bool skip = tok[0] == "class" || tok[0][0] == '#' ;

      bool include = true;
      bool invalid = false;
      
      if ( ! ( header || skip ) )
	{

	  // skip '.' rows in annots
	  if ( tok[ col1-1 ] == "." || tok[ col2-1 ] == "." )
	    invalid = true;
	  
	  double t1 = atof( tok[ col1-1 ].c_str() );
	  double t2 = atof( tok[ col2-1 ].c_str() );
	  
	  if ( t1 > t2 )
	    std::cerr << " *** bad region: " << line << "\n";
	  else
	    {	      
	      reg_t r( t1 , t2 );	      
	      include = region_is_included( r , regions );
	    }
	  
	}
      
      //
      // print 
      //

      
      if ( header || skip )
	{
	  std::cout << line;
	  if ( append_col && header ) std::cout << "\t" << newcol;
	  std::cout << "\n";
	}
      else if ( ! invalid ) 
	{
	  ++data_line;
	  if ( include ) ++included_data_line;
	  
	  if ( append_col ) 
	    std::cout << line << "\t" << (int)include << "\n";
	  else
	    {
	      if ( exclude_matches && ! include )
		std::cout << line << "\n";
	      else if ( include_matches && include )
	 	std::cout << line << "\n";
	    }
	  
	}
      
    }
  
  std::cerr << included_data_line << " matched, " << data_line -included_data_line << " unmatched, of "
	    << data_line << " valid lines in total\n";
  
  std::exit(0);
}


bool region_is_included( const reg_t & r , const std::set<reg_t> & d )
{

  // does region 'r' overlap any region in 'd'?
  // n.b. 'd' will be flattened
  if ( d.size() == 0 ) return false;
  
  // find the first annot not before (at or after) this region
  std::set<reg_t>::const_iterator lwr = d.lower_bound( r );

  // in case starts at the same position
  if ( lwr != d.end() )
    {      
      if ( lwr->start <= r.stop && lwr->stop >= r.start ) return true;
    }
  
  // slide back one (i.e. as region starting before 'r' could still
  // overlap 'r'
  
  if ( lwr != d.begin() )
    {
      --lwr;
      if ( lwr->start <= r.stop && lwr->stop >= r.start ) return true;
    }
  
  // because we flattened list, we do not have to go back any further (i.e.
  // end of prior will be before start of this one)
    
  return false;
}


std::set<reg_t> flatten( const std::set<reg_t> & x )
{

  std::set<reg_t> m;
  
  if ( x.size() == 0 ) return m;
  if ( x.size() == 1 ) return x;
  
  reg_t curr = *x.begin();
  
  std::set<reg_t>::const_iterator xx = x.begin();
  while ( xx != x.end() )
    {
      // because of +1 end encoding, contiguous regions will (a,b) ... (b,c) 
      const reg_t & pro = *xx;
      if ( pro.start > curr.stop )
	{	      
	  m.insert( curr );	      
	  	      
	  // and update the current
	  curr = pro;
	}
      else // expand background as needed
	{
	  if ( pro.stop > curr.stop ) curr.stop = pro.stop;
	}

      // consider next element
      ++xx;
    }
  
  // add final element
  m.insert( curr );

  return m;
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
