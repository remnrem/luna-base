
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
#include <fstream>
#include <cstdlib>
#include <set>


struct obs_t { 
  obs_t() { } 
  obs_t( const std::string & obs , const std::string & indiv ) : obs(obs) , indiv(indiv) { } 
  std::string indiv;
  std::string obs;
  bool operator<( const obs_t & rhs ) const { if ( indiv == rhs.indiv ) return obs < rhs.obs; return indiv < rhs.indiv; } 
};

#include "helper/helper.h"
#include "miscmath/crandom.h"

int main(int argc, char ** argv)
{
  if ( argc != 3 ) Helper::halt("tspicker list N < input" );
  int n;
  if ( ! Helper::str2int( argv[2] , &n ) ) Helper::halt("tspicker list N < input" );
  std::string filename = argv[1];
  if ( ! Helper::fileExists( filename ) ) Helper::halt("tspicker list N < input" );
  
  std::set<obs_t> obs;
  std::vector<obs_t> nobs;
  
  std::ifstream IN1( filename.c_str() , std::ios::in );
  while ( ! IN1.eof() ) 
    {
      obs_t o;
      IN1 >> o.obs >> o.indiv;
      if ( IN1.eof() ) continue;
      if ( o.obs == "" ) continue;
      obs.insert( o );
      nobs.push_back( o );
    }
  IN1.close();
  
  std::cerr << "read " << obs.size() << " obs\n";

  if ( obs.size() != nobs.size() ) Helper::halt( "duplicates found" );

  const int s = obs.size();
  
  // select N from these 
  CRandom::srand(time(0));

  // clear this, to use as final N-selected list
  obs.clear();

  std::set<int> incl;
  
  if ( s < n ) Helper::halt( "not enough obs" );

  int cnt = 0;
  while ( cnt < n )
    {
      int r = CRandom::rand( s );
      // already in?
      if ( incl.find( r ) != incl.end() ) continue;
      
      //      std::cout << nobs[r].obs << "\t" << nobs[r].indiv << "\n";
      
      incl.insert( r );
      obs.insert( nobs[r] );
      ++cnt;
      
    }
  
  std::cerr << "selecting " << n << " from " << s << " obs\n";
  
  //
  // Now read through STDIN, only selecting lines that match 
  //
  
  while ( ! std::cin.eof() )
    {
      std::string line;
      Helper::safe_getline( std::cin , line );
      if ( line == "" ) continue;
      if ( std::cin.eof() ) break;
      std::vector<std::string> tok = Helper::parse( line , "\t" );
      obs_t o( tok[0] , tok[1] );
      if ( obs.find(o) != obs.end() ) std::cout << line << "\n";
    }
  
  
  std::exit(0);
}

