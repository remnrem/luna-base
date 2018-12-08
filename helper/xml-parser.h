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
