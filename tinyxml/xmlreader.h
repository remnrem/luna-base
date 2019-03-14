
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

#ifndef __XMLREADER_H__
#define __XMLREADER_H__

#include "tinyxml.h"
#include "helper/helper.h"
#include <string>
#include <vector>
#include <map>


struct element_t;
struct attr_t;

class XML
{

  friend struct element_t;

 public:
  
 XML( const std::string & f )
   : doc( f ) 
  {
    filename = f;          
    is_valid = doc.LoadFile();      
    if ( is_valid ) parse( &doc );
  }
  
  ~XML();

  bool valid() const { return is_valid; }  
  void dump();
  void dumper( element_t * e );
  void dump_history( element_t * e , std::vector<std::string> * );
  
  std::vector<element_t*> children( const std::string & e );

 private:  
  std::string filename;  
  element_t * root;
  TiXmlDocument doc;  
  bool is_valid;  
  attr_t parse_attr( TiXmlElement* pElement );  
  void parse( TiXmlNode * pParent , element_t * pGrandParent = NULL );

  static void finder( element_t * e , const std::string & n , std::vector<element_t*> ** r );
  
};

struct attr_t 
{
  void add( const std::string & key , const std::string & value )
  {
    std::pair<std::string,std::string> p( key , value ) ;
    alist.push_back( p );
    amap[ key ] = value;
  }
  std::vector<std::pair<std::string,std::string> > alist;
  std::map<std::string,std::string> amap;
  int size() const { return alist.size(); }
  std::string key(const int i) const { return alist[i].first; }
  std::string value(const int i) const { return alist[i].second; }
  std::string value(const std::string & key ) { return amap[key]; }
};
  
struct element_t
{
  
  friend class XML;

  element_t( const std::string & name , element_t * parent = NULL )
  : parent(parent) , name(name) , value("") 
  {
    if ( parent ) parent->child.push_back( this );
  }
  
  ~element_t( ) 
  {
    for (int c=0;c<child.size();c++) delete child[c];    
  }
  
  element_t * parent;

  std::vector<element_t*> child;

  element_t * operator()( const std::string & n ) const
  {
    for ( int c=0;c<child.size();c++) if ( Helper::iequals( child[c]->name , n ) ) return child[c];
    return NULL;
  }
  
  std::vector<element_t*> children( const std::string & n ) 
  {
    std::vector<element_t*> * r = NULL;
    XML::finder( this , n , &r );  
    if ( r == NULL ) { std::vector<element_t*> dummy; return dummy; } 
    return *r;
  }
  
  // <name attr1=1 attr2=text>Value is here</name>
  std::string name;
  std::string value;
  attr_t attr;
};

#endif
