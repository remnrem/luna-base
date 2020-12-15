

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

#ifndef __CACHE_H__
#define __CACHE_H__

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include "helper/helper.h"

// a temporary store that can be used for different commands to communicate to each other
//  SPINDLES ... cache=label1
//  SPINDLES ... cache=label2
//  SO .. cache=label3
//  SO .. cache=label4
//  COUPL spindles=label1 so=label3
//  COUPL spindles=label1 so=label4
//  etc.

// where the different commands know what to store/expect
//  e.g. SPINDLES saves spindle peak time-points and all wavelet coefficients
//       SO saves SO time-points and slow phase
//  the cache is also defined by stratifiers, e.g. label1 / F=11 / CH == C3
//   





void ctest();

struct ckey_t {
  
  ckey_t( const std::string & name , const std::map<std::string,std::string> & stratum )
    : name( name ) , stratum( stratum )
  {
  }
  
  ckey_t( const std::string & name ) : name( name )
  {
  }
  
  void add( const std::string & key , const std::string & val )
  {
    stratum[key] = val;
  }
  
  void add( const std::string & key , const int & val )
  {
    stratum[key] = Helper::int2str( val );
  }
  
  void add( const std::string & key , const double & val )
  {
    stratum[key] = Helper::dbl2str( val );
  }
  
  void add( const std::string & key , const bool val )
  {
    stratum[key] = Helper::int2str( val );
  } 
  
  std::string name;

  std::map<std::string,std::string> stratum;
  
  bool operator<(const ckey_t & rhs ) const {
    if ( name < rhs.name ) return true;
    if ( name > rhs.name ) return false;
    if ( stratum.size() < rhs.stratum.size() ) return true;
    if ( stratum.size() > rhs.stratum.size() ) return false;
    std::map<std::string,std::string>::const_iterator ii = stratum.begin();
    std::map<std::string,std::string>::const_iterator jj = rhs.stratum.begin();	
    while ( ii != stratum.end() )
      {
	if ( ii->first < jj->first ) return true;
	if ( ii->first > jj->first ) return false;
	if ( ii->second < jj->second ) return true;
	if ( ii->second > jj->second ) return false;
	++ii;
	++jj;
      }
    return false;
  }
  
};
  
  
template<class T>
struct cache_t {

  //  cache_t() { } 

  cache_t( const std::string & name ) : name(name) { }   

  std::string name;

  std::map<ckey_t,std::vector<T> > store;

  // member functions
  
  void add( const ckey_t & key , const std::vector<T> & value )
  {
    store[key] = value;
  }

  void add( const ckey_t & key , const T & value )
  {
    std::vector<T> t(1);
    t[0] = value;
    add( key , t );
  }

  void clear()
  {
    store.clear();
  }

  // get all keys matching a particular label
  std::set<ckey_t> keys( const std::string & n ) const {
    std::set<ckey_t> k;
    typename std::map<ckey_t,std::vector<T> >::const_iterator ii = store.begin();
    while ( ii != store.end() )
      {
	if ( ii->first.name == n ) k.insert( ii->first );
	++ii;
      }
    return k;
  }

  std::vector<T> size( const ckey_t & key ) const {
    typename std::map<ckey_t,std::vector<T> >::const_iterator ii = store.find( key );
    if ( ii == store.end() ) return 0;
    return ii->second.size();
  }
  
  
  std::vector<T> fetch( const ckey_t & key ) const {
    typename std::map<ckey_t,std::vector<T> >::const_iterator ii = store.find( key );
    if ( ii == store.end() )
      {
	std::vector<T> dummy; return dummy; 
      }
    return ii->second;
  }

  
  std::string print() const {
    std::stringstream oo;
    oo << "cache: " << name << "\n";
    typename std::map<ckey_t,std::vector<T> >::const_iterator ss = store.begin();
    while ( ss != store.end() )
      {
	oo << "\t" << ss->first.name << "\n";
	std::map<std::string,std::string>::const_iterator kk = ss->first.stratum.begin();
	while ( kk != ss->first.stratum.end() )
	  {
	    oo << "\t" << kk->first << " --> " << kk->second << "\n";
	    ++kk;
	  }
	oo << "\tdata: " << ss->second.size() << " element vector\n";
	++ss;
      }
    return oo.str();
  }
  
  
};



//
// hihgest level interface for caches (member of timeline_t)
//

struct caches_t {

  std::map<std::string, cache_t<int> > cache_int;
  std::map<std::string, cache_t<double> > cache_num;  
  std::map<std::string, cache_t<std::string> > cache_str;  
  std::map<std::string, cache_t<uint64_t> > cache_tp;

  void clear() {
    cache_int.clear();
    cache_num.clear();
    cache_str.clear();
    cache_tp.clear();
  }
  
  bool has_int( const std::string & n ) const { return cache_int.find(n) != cache_int.end(); }
  bool has_str( const std::string & n ) const { return cache_str.find(n) != cache_str.end(); }
  bool has_num( const std::string & n ) const { return cache_num.find(n) != cache_num.end(); }
  bool has_tp( const std::string & n ) const { return cache_tp.find(n) != cache_tp.end(); } 
  

  //
  // fetch, and create if does not exist
  //

  cache_t<int> * find_int( const std::string & n )
  {
    if ( ! has_int( n ) ) cache_int.insert( std::pair<std::string,cache_t<int> >( n , cache_t<int>( n ) ) );;
    return &(cache_int.find( n )->second );
  }

  cache_t<std::string> * find_str( const std::string & n )
  {
    if ( ! has_str( n ) ) cache_str.insert( std::pair<std::string,cache_t<std::string> >( n , cache_t<std::string>( n ) ) );;
    return &(cache_str.find( n )->second );
  }

  cache_t<double> * find_num( const std::string & n )
  {
    if ( ! has_num( n ) ) cache_num.insert( std::pair<std::string,cache_t<double> >( n , cache_t<double>( n ) ) );;
    return &(cache_num.find( n )->second );
  }

  cache_t<uint64_t> * find_tp( const std::string & n )
  {
    if ( ! has_tp( n ) ) cache_tp.insert( std::pair<std::string,cache_t<uint64_t> >( n , cache_t<uint64_t>( n ) ) );;
    return &(cache_tp.find( n )->second );
  }


  
};


#endif



