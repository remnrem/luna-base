
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

#include "xmlreader.h"


std::vector<element_t*> XML::children( const std::string & n )
{
  // search entire documenet to first instance of this type
  // then return all kids
  std::vector<element_t*> * r = NULL;
  finder( root , n , &r );  
  if ( r == NULL ) { std::vector<element_t*> dummy; return dummy; } 
  return *r;
}

void XML::finder( element_t * e , const std::string & n , std::vector<element_t*> ** r ) 
{
  // inefficient, but assume we'll usually be looking for label near the start of 
  // the file (i.e. that earlier siblings won't have lots of children
  if ( Helper::iequals( e->name , n ) ) *r = &(e->child);
  for (int c=0;c<e->child.size();c++) finder( e->child[c] , n , r );  
}


void XML::dump()
{
  if ( root ) dumper(root);
}


void XML::dump_history( element_t * e , std::vector<std::string> * history )
{
  if ( e->parent == NULL ) return;
  history->push_back( e->parent->name );
  dump_history( e->parent , history );
  return;
}

void XML::dumper( element_t * e )
{  
  // history
  std::vector<std::string> history;
  dump_history( e , &history);
  std::vector<std::string>::reverse_iterator hh = history.rbegin();
  while ( hh != history.rend() ) { std::cout << *hh << "|" ; ++hh; } 
  std::cout << e->name << " = " << e->value << "\t[ ";
  int na = e->attr.size();
  for (int i=0; i<na;i++) std::cout << e->attr.key(i) << "=" << e->attr.value(i) << " ";
  std::cout << "]\n";
  // kids
  for (int c=0;c<e->child.size();c++) dumper( e->child[c] );  
}

attr_t XML::parse_attr( TiXmlElement* pElement )
{
  attr_t a;  
  if ( pElement == NULL) return a;    
  TiXmlAttribute* pAttrib=pElement->FirstAttribute();    
  while ( pAttrib )
    {
      a.add( pAttrib->Name() , pAttrib->Value() );
      pAttrib=pAttrib->Next();
    }
  return a;
}


void XML::parse( TiXmlNode* pNode , element_t * pGrandParent )
{

  if ( pNode == NULL ) return;

  element_t * current = pGrandParent;  
  TiXmlNode * pChild;  
  TiXmlText * pText;
  
  switch ( pNode->Type() )
    {
      
      //
      // root node
      //
      
    case TiXmlNode::TINYXML_DOCUMENT:
      root = new element_t( "Document" );
      current = root;
      break;
      
      //
      // element
      //
      
    case TiXmlNode::TINYXML_ELEMENT:      
      current = new element_t( pNode->Value() , pGrandParent );     
      // parse attributes
      current->attr = parse_attr( pNode->ToElement() );      
      break;

      
      //
      // Text/values
      //
      
    case TiXmlNode::TINYXML_TEXT:
      if ( pGrandParent ) 
	pGrandParent->value = pNode->ToText()->Value();     	
      break;

      //
      // Ignore everything else
      //

    case TiXmlNode::TINYXML_COMMENT:
      //printf( "Comment: [%s]", pNode->Value());
      break;

    case TiXmlNode::TINYXML_UNKNOWN:
      //printf( "Unknown" );
      break;

    case TiXmlNode::TINYXML_DECLARATION:
      //printf( "Declaration" );
      break;

    default:
      break;
    }


  //
  // add all children
  //

  if ( current ) 
    {
      for ( pChild = pNode->FirstChild(); 
	    pChild != 0; 
	    pChild = pChild->NextSibling()) 
	parse( pChild, current );	
    }
}

