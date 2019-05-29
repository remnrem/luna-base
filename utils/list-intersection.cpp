
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

#include "intervals/intervals.cpp"
#include "helper/helper.h"

#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <cstring>
#include <fstream>

extern globals global;

int main(int argc , char ** argv )
{
  
  global.init_defs();

  global.api();

  if ( argc < 3 ) Helper::halt( "usage: intersect list1 list2 -b bglist -x exclude-list -w window -t threshold (0..1) {-s} " );
  
  std::string list1, list2;

  bool has_bg = false, xlist = false; 
  bool sec_mode = false;
  
  std::string bglist = "";
  
  // param
  double threshold = 0.5;  // -t
  double window = 0;       // -w
  int listed = 0;

  for (int i=1;i<argc;i++)
    {
      
      if ( strcmp( argv[i] , "-t" ) == 0 ) 
	{
	  if ( i < argc-1 && Helper::str2dbl( argv[i+1] , &threshold) ) { ++i; continue; }
	  else Helper::halt("bad -t option"); 
	}
      
      else if ( strcmp( argv[i] , "-w" ) == 0 ) 
	{
	  if ( i < argc-1 && Helper::str2dbl( argv[i+1] , &window) ) { ++i; continue; }
	  else Helper::halt("bad -w option"); 
	}
      
      else if ( strcmp( argv[i] , "-b" ) == 0 ) 
	{
	  has_bg = true; xlist = false; 
	  bglist = argv[i+1]; ++i; continue;
	}

     else if ( strcmp( argv[i] , "-x" ) == 0 ) 
	{
	  has_bg = true; xlist = true;
	  bglist = argv[i+1]; ++i; continue;	  
	}
     else if ( strcmp( argv[i] , "-s" ) == 0 ) 
       {
	 // expect doubles (seconds) and not TP-units
	 sec_mode = true;
       }
      
      else 
	{
	  if ( listed > 2 ) Helper::halt( "too many lists" );
	  ++listed;
	  if ( listed == 1 ) list1 = argv[i];
	  if ( listed == 2 ) list2 = argv[i];
	}
    }

  uint64_t window_msec = window <= 0 ? 0 : window * globals::tp_1sec;

  // lists

  std::set<interval_t> bg, l1, l2;
  std::map<interval_t,std::vector<std::string> > lines1, lines2;

  
  // background intervals?

  if ( has_bg )
    {

      if ( ! Helper::fileExists( bglist ) ) Helper::halt( "could not find " + bglist );

      std::cerr << "opening bg-list " << bglist << "\n";
      std::ifstream IN1( bglist.c_str() , std::ios::in );
      while ( !IN1.eof() ) 
	{      
	  std::string line;      
	  Helper::safe_getline( IN1 , line );      
	  if ( IN1.eof() ) continue;
	  std::vector<std::string> tok = Helper::parse( line , "\t" );
	  const int n = tok.size();
	  
	  if ( n < 2 ) Helper::halt( "bad line in bg-list:\n" + line );
	  

	  uint64_t a1, b1;

	  if ( sec_mode )
	    {
	      double a, b;
	      if ( ! Helper::str2dbl( tok[0] , &a ) ) Helper::halt( "bad format" );
	      if ( ! Helper::str2dbl( tok[1] , &b ) ) Helper::halt( "bad format" );
	      
	      a1 = a * globals::tp_1sec;
	      b1 = b * globals::tp_1sec;
	    }
	  else
	    {
	      if ( ! Helper::str2int64( tok[0] , &a1 ) ) Helper::halt( "bad format" );
	      if ( ! Helper::str2int64( tok[1] , &b1 ) ) Helper::halt( "bad format" );
	    }
	  
	  
	  interval_t i( a1, b1 );
	  
	  bg.insert( i );
	}
      IN1.close();
      std::cerr << "read " << bg.size() << " background elements to " 
		<< ( xlist ? "exclude" : "include" ) << "\n";
    }


  // read lists: expecting FTR format -- will save rest of line, just take first two fields

  if ( ! Helper::fileExists( list1 ) ) Helper::halt( "could not find " + list1 );

  std::ifstream IN1( list1.c_str() , std::ios::in );
  while ( !IN1.eof() ) 
    {      
      std::string line;      
      Helper::safe_getline( IN1 , line );      
      if ( IN1.eof() ) continue;
      std::vector<std::string> tok = Helper::parse( line , "\t" );
      const int n = tok.size();
      if ( n < 2 ) Helper::halt( "bad line in list 1:\n" + line );
      
      uint64_t a1, b1;
      
      if ( sec_mode )
	{
	  double a,b;
	  
	  if ( ! Helper::str2dbl( tok[0] , &a ) ) Helper::halt( "bad format" );
	  if ( ! Helper::str2dbl( tok[1] , &b ) ) Helper::halt( "bad format" );
	  
	  a1 = a * globals::tp_1sec;
	  b1 = b * globals::tp_1sec;
	  
	}
      else
	{
	  if ( ! Helper::str2int64( tok[0] , &a1 ) ) Helper::halt( "bad format" );
	  if ( ! Helper::str2int64( tok[1] , &b1 ) ) Helper::halt( "bad format" );
	}
      
      interval_t i( a1, b1 );
      
      //save line
      l1.insert( i );
      lines1[ i ] = tok;

    }
  IN1.close();


  if ( ! Helper::fileExists( list2 ) ) Helper::halt( "could not find " + list2 );

  std::ifstream IN2( list2.c_str() , std::ios::in );
  while ( !IN2.eof() ) 
    {      
      std::string line;      
      Helper::safe_getline( IN2 , line );      
      if ( IN2.eof() ) continue;
      std::vector<std::string> tok = Helper::parse( line , "\t" );
      const int n = tok.size();
      if ( n < 2 ) Helper::halt( "bad line in list 2:\n" + line );

      uint64_t a1, b1;

      if ( sec_mode )
	{
	  double a,b;
	  
	  if ( ! Helper::str2dbl( tok[0] , &a ) ) Helper::halt( "bad format" );
	  if ( ! Helper::str2dbl( tok[1] , &b ) ) Helper::halt( "bad format" );
	  
	  a1 = a * globals::tp_1sec;
	  b1 = b * globals::tp_1sec;
	  
	}
      else
	{
	  if ( ! Helper::str2int64( tok[0] , &a1 ) ) Helper::halt( "bad format" );
	  if ( ! Helper::str2int64( tok[1] , &b1 ) ) Helper::halt( "bad format" );
	}
      
      interval_t i( a1, b1 );
      
      //save line
      l2.insert( i );
      lines2[ i ] = tok;

    }
  IN2.close();


  //  std::cerr << "lists contained " << l1.size() << " and " << l2.size() << " intervals\n";


  // need to first intersect with the background?
  
  if ( has_bg )
    {
      int l1sz = l1.size(); int l2sz = l2.size();
      
      std::set<interval_t> ba, bb, cons, uns, oa, ob;
      std::set<interval_t> n1, n2;
 
      // prune list 1 (0,0 = any overlap, 0 window  around intervals)
      int olap = interval_t::intersect( l1, bg , &ba, &bb, &cons, &uns, &oa, &ob , 0 , 0 );
      
      std::set<interval_t>::const_iterator ii1 = l1.begin();
      while ( ii1 != l1.end() )
	{
	  // overlap?
	  bool overlaps = ba.find( *ii1 ) != ba.end();
	  
	  if ( xlist && ! overlaps )  n1.insert( *ii1 );
	  else if ( overlaps && ! xlist )  n1.insert( *ii1 );
	  
	  ++ii1;
	}
      
      ba.clear(); bb.clear(); cons.clear(); uns.clear(); oa.clear(); ob.clear();

      // prune list 1 (0,0 = any overlap, 0 window  around intervals)
      olap = interval_t::intersect( l2, bg , &ba, &bb, &cons, &uns, &oa, &ob , 0 , 0 );
      
      std::set<interval_t>::const_iterator ii2 = l2.begin();
      while ( ii2 != l2.end() )
	{
	  // overlap?
	  bool overlaps = ba.find( *ii2 ) != ba.end();
	  
	  if ( xlist && ! overlaps )  n2.insert( *ii2 );
	  else if ( overlaps && ! xlist )  n2.insert( *ii2 );
	  
	  ++ii2;
	}
           
      // copy over
      l1 = n1;
      l2 = n2;
      
      std::cerr << "applying background list:\n" << list1 << " retains " << l1.size() << " of " << l1sz << " intervals\n";
      std::cerr << list2 << " retains " << l2.size() << " of " << l2sz << " intervals\n";
    } 
  
//    std::cerr << "read " << l1.size() << " interals from " << list1 << "\n";
//    std::cerr << "read " << l2.size() << " interals from " << list2 << "\n";
  
  std::set<interval_t> ba, bb, cons, uns, oa, ob;
  
  // intersect
  
  int olap = interval_t::intersect( l1, l2 , &ba, &bb, &cons, &uns, &oa, &ob , threshold , window_msec );
  
  double olapa = ba.size() / (double)l1.size() ;
  double olapb = bb.size() / (double)l2.size() ;

  std::cerr << "# intervals : " << l1.size() << "\t" << l2.size() << "\n";
  std::cerr << "# overlap   : " << ba.size() << "\t" << bb.size() << "\n";
  
  std::cerr << "p(overlap)  : ";
  if ( l1.size() > 0 ) std::cerr << olapa ; else std::cerr << "n/a";
  std::cerr << "\t";
  if ( l2.size() > 0 ) std::cerr << olapb ; else std::cerr << "n/a";
  std::cerr << "\n";
  
  if ( l1.size() > 0 && l2.size() > 0 )
    std::cerr << "average p(overlap) : " << ( olapa + olapb ) / 2.0 << "\n";
  
  //
  // stats
  //

  // list overlapping spindles 
  
  //   bool olap = oa.find( ii->tp ) == oa.end();
  //   bool olap = ob.find( ii->tp ) == ob.end();
  
  std::set<interval_t>::const_iterator ii1 = l1.begin();
  while ( ii1 != l1.end() )
    {
      const interval_t & interval = *ii1;
      
      // overlap?
      bool overlaps = ba.find( interval ) != ba.end();
      
      std::cout << "A\t"
		<< overlaps;      
      std::vector<std::string> & lines = lines1[ interval ];
      for (int i=0;i<lines.size();i++) 
	{
	  std::cout << "\t" << lines[i];
	  std::vector<std::string> tok2 = Helper::parse( lines[i] , "=" );
	  //if ( tok2.size() == 2 ) std::cout << "\t" << tok2[1];
	}
      std::cout << "\n";
      
      ++ii1;
    }
  
  
  //  B list

  std::set<interval_t>::const_iterator ii2 = l2.begin();
  while ( ii2 != l2.end() )
    {
      const interval_t & interval = *ii2;
      
      // overlap?
      bool overlaps = bb.find( interval ) != bb.end();
      
      std::cout << "B\t"
		<< overlaps;      
      std::vector<std::string> & lines = lines2[ interval ];
      for (int i=0;i<lines.size();i++) 
	{
	  std::cout << "\t" << lines[i];
	  std::vector<std::string> tok2 = Helper::parse( lines[i] , "=" );
	  //if ( tok2.size() == 2 ) std::cout << "\t" << tok2[1];
	}
      std::cout << "\n";
      
      ++ii2;
    }
  
  std::exit(0);
}
