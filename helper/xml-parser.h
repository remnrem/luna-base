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

#ifndef __XML_PARSER_H__
#define __XML_PARSER_H__

#include "tinyxml/tinyxml.h"

// wrapper around TinyXML

class MyXML
{
  
 public:
  
  MyXML( const std::string & filename ) : doc(filename) 
  { 
    //
  } 
    
  bool load() { return doc.LoadFile(); } 
  
  void dump() { dump_to_stdout( &doc ); }
  
 private:
  
  TiXmlDocument doc;

  static const unsigned int NUM_INDENTS_PER_SPACE = 2;
  
  const char * getIndent( unsigned int numIndents );
  const char * getIndentAlt( unsigned int numIndents );
  int dump_attribs_to_stdout(TiXmlElement* pElement, unsigned int indent);
  void dump_to_stdout( TiXmlNode* pParent, unsigned int indent = 0 );
  
};

#endif
