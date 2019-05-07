
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


#include "sstore/sstore.h"
#include <iostream>
#include "helper/helper.h"
#include "helper/logger.h"
#include <vector>
#include <set>
#include <map>

extern logger_t logger;

int main(int argc , char ** argv )
{
  
  logger.off();

  std::map<std::string,int> lvls;
  for (int i=1;i<argc;i++)
    lvls[ argv[i] ] = -1 ;
    
  // get header line
  
  std::string line;
  Helper::safe_getline( std::cin , line  );
  std::vector<std::string> hdr = Helper::parse( line , "\t" );
  const int n = hdr.size();
  int eidx = -1;
  int iidx = -1;
  int chidx = -1;
  int startidx = -1;
  int stopidx = -1;
  bool has_levels = false;

  std::map<std::string,int> vars;

  for (int i=0;i<n;i++)
    {
      if      ( hdr[i] == "E" ) eidx = i;
      else if ( hdr[i] == "START" ) startidx = i;
      else if ( hdr[i] == "STOP" ) stopidx = i;
      else if ( hdr[i] == "ID" ) iidx = i;
      else if ( hdr[i] == "CH" ) chidx = i;
      else if ( hdr[i] == "CHS" ) chidx = i;       // use this as an alternative for COH analyses
      else if ( lvls.find( hdr[i] ) != lvls.end() ) { lvls[ hdr[i] ] = i; has_levels = true; } 
      else vars[ hdr[i] ] = i;
    }

  bool epochs = eidx != -1;
  bool intervals = startidx != -1 & stopidx != -1;
  if ( iidx == -1 ) Helper::halt( "no ID column" );
  if ( epochs && intervals ) Helper::halt( "cannot have both intervals and epochs" );
  
  bool base = ! ( epochs || intervals );

  bool has_channels = chidx != -1;
  
  //
  // read and process each non-header row
  //

  while ( ! std::cin.eof() )
    {
      
      std::string line;
      Helper::safe_getline( std::cin , line );
      if ( std::cin.eof() ) break;
      if ( line == "" ) break;
      
      std::vector<std::string> tok = Helper::parse( line , "\t" );
      
      // level?
      std::string lvl = ".";
      if ( has_levels ) 
	{
	  lvl = "";
	  std::map<std::string,int>::const_iterator ii = lvls.begin();
	  while ( ii != lvls.end() ) 
	    {
	      if ( ii != lvls.begin() ) lvl += ";";
	      lvl += ii->first + "=" + tok[ ii->second ];
	      
	      ++ii;	      
	    }	  
	}
      
      // channels?
      std::string ch = has_channels ? tok[ chidx ] : ".";
      
      std::map<std::string,int>::const_iterator ii = vars.begin();
      while ( ii != vars.end() )
	{
	  
	  double d;
	  bool is_double = Helper::str2dbl( tok[ ii->second ] , &d );
	  
	  std::cout << ii->first << "\t"
		    << lvl << "\t"
		    << ch << "\t" ;
	  
	  if ( epochs ) 
	    {
	      std::cout << tok[ eidx] << "\t";
	    }
	  else if ( intervals )
	    {
	      std::cout << tok[ startidx ] << "\t" 
			<< tok[ stopidx ] << "\t";
	    }
	  
	  std::cout << ( is_double ? 1 : 0 ) << "\t"
		    << tok[ ii->second ] << "\n";
	  
	  ++ii;
	}

    }

  std::exit(0);
}
