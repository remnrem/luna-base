
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

#include "timeline/cache.h"

#include "db/db.h"
#include "helper/logger.h"

extern writer_t writer;
extern logger_t logger;

void ctest()
{
  
  writer.level( "L1" , "F1" );
  
  writer.level( 123 , "FFE" );
  writer.epoch( 222 );
  
  cache_t<double> cache( "my1" );

  ckey_t ckey1( "y" , writer.faclvl() ) ;
  ckey_t ckey2( "z" , writer.faclvl() ) ;

  std::vector<double> y( 10 , 22 );
  std::vector<double> z( 10 , 23 );
  
  cache.add( ckey1 , y );
  cache.add( ckey2 , z );
  
  writer.unlevel();

  std::cout << cache.print();
  
}


void caches_t::load( const std::string & filename )
{
  
  std::ifstream IN1( filename.c_str() , std::ios::in );
  

  // 
  // format : 
  //
  // cache: peaks[int]
  // strata: fac=lvl
  // strata: clear
  // value: var=123.456
  

  std::map<std::string,std::string> curr_strata;
  
  cache_t<std::string> * str_cache = NULL ;
  cache_t<double> * num_cache = NULL ;
  cache_t<int> * int_cache = NULL ;
  cache_t<uint64_t> * tp_cache = NULL ;
  
  int cnt = 0 ;

  while ( ! IN1.eof() )
    {
      
      std::string line;

      Helper::safe_getline( IN1 , line );
      
      if ( IN1.eof() ) break;
      if ( line == "" ) continue;
      
      std::vector<std::string> tok = Helper::parse( line , "\t " );
      if ( tok.size() != 2 ) Helper::halt( "problem with cache format: " + line );
      
      if ( tok[0] == "cache:" ) 
	{
	  std::vector<std::string> tok2 = Helper::parse( line , "[]" );

	  str_cache = NULL;
	  int_cache = NULL;
	  num_cache = NULL;
	  tp_cache = NULL;

	  if ( tok2.size() != 2 ) Helper::halt( "problem with cache format: " + line );
	  if ( tok2[1] == "int" )
	    int_cache = find_int( tok2[0] );
	  else if ( tok2[1] == "num" )
	    num_cache = find_num( tok2[0] );
	  else if ( tok2[1] == "str" )
	    str_cache = find_str( tok2[0] );
	  else if ( tok2[1] == "tp" )
	    tp_cache = find_tp( tok2[0] );
	  else
	    Helper::halt( "problem with cache format: " + line );
	  logger << "reading into " << tok2[0] << "\n";
	}
      else if ( tok[0] == "strata:" )
	{
	  if ( tok[1] == "clear" ) curr_strata.clear();
	  else {
	    std::vector<std::string> tok2 = Helper::parse( tok[1], "=" );
	    if ( tok2.size() != 2 )  Helper::halt( "problem with cache format: " + line);
	    curr_strata[ tok2[0] ] = tok2[1];
	  }
	}
      else if ( tok[0] == "value:" )
	{
	  std::vector<std::string> tok2 = Helper::parse( tok[1], "=" );
	  if ( tok2.size() !=2 )  Helper::halt( "problem with cache format: " + line);
	  
	  if ( num_cache )
	    {
	      double num_value;
	      if ( ! Helper::str2dbl( tok2[1] , &num_value ) ) Helper::halt( "problem with cache format: " + line);
	      num_cache->add( ckey_t( tok2[0] , curr_strata ) , num_value ) ;
	      std::cout << " adding " << tok2[0] << " --> " << tok2[1] << "\n";
	    }
	  else if ( int_cache )
	    {
	      int int_value;
	      if ( ! Helper::str2int( tok2[1] , &int_value ) ) Helper::halt( "problem with cache format: " + line);
	      int_cache->add( ckey_t( tok2[0] , curr_strata ) , int_value ) ;
	    }
	  else if ( str_cache )
	    {
	      str_cache->add( ckey_t( tok2[0] , curr_strata ) , tok2[1] );
	    }
	  else if ( tp_cache )
	    {
	      uint64_t tp_value;
	      if ( ! Helper::str2int64( tok2[1] , &tp_value ) ) Helper::halt( "problem with cache format: " + line);
	      tp_cache->add( ckey_t( tok2[0] , curr_strata ) , tp_value ) ; 
	    }
	  else
	    Helper::halt( "problem with cache format: " + line);
	  
	  ++cnt;

	}
      else
	Helper::halt( "problem with cache format: " + line);

      // next row
    }
  
  IN1.close();

  logger << "  read " << cnt << " values from " << filename << "\n";

  std::cout << " print \n\n" << num_cache->print() << "\n\n---\n";
  
}

