
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

#include "utils/cgi-utils.h"
#include "luna.h"

extern writer_t writer;

extern globals global;

void input_page( const std::map<std::string,std::string> & );
void output_page( const std::map<std::string,std::string> & );
  
int main( int argc , char ** argv )
{

  // set up Luna 
  global.init_defs();
  global.api();

  // clear all existing mappings
  nsrr_t::clear();
  
  // set up HTML output
  html_write_headers( "NAP mapper" );
  
  std::cout << "<h2 style=\"color:navy;font-family:georgia\">NAP channel/annotation mapper</h2>";
  std::cout << "mappings: <a href=\"https://gitlab-scm.partners.org/zzz-public/nsrr/-/tree/master/common/resources\">gitlab repo</a>";  
  std::cout << "<hr>";
  
  // get arguments
  std::map<std::string,std::string> vars = fetch_cgi();

  // inputs?
  if ( vars.size() == 0 || vars.find( "inp" ) != vars.end() )
    {
      input_page( vars );
    }
  else
    {
      output_page( vars );
    }
       
  html_write_footer();

}


void input_page( const std::map<std::string,std::string> & vars )
{

  std::string str1, str2, str3, str4;

  if ( vars.find( "f1" ) != vars.end() ) str1 = Helper::search_replace( vars.find( "f1" )->second , '|' , '\n' );
  if ( vars.find( "f2" ) != vars.end() ) str2 = Helper::search_replace( vars.find( "f2" )->second , '|' , '\n' );
  if ( vars.find( "f3" ) != vars.end() ) str3 = Helper::search_replace( Helper::search_replace( vars.find( "f3" )->second , '|' , '\n' ) , '=' , ' ' );
  if ( vars.find( "f4" ) != vars.end() ) str4 = Helper::search_replace( Helper::search_replace( vars.find( "f4" )->second , '|' , '\n' ) , '=' , ' ');
        
  std::cout << "<form action=\"/cgi-bin/cgi-mapper\" method=\"post\">"    
	    << "<table width=\"100%\" border=0><tr width=\50%\"><td>"
	    << "<label for=\"f1\">Channels</label><br>"
	    << "<textarea style=\"width: 100%; max-width: 100%; resize: none; font-family: Courier New\" id=\"f1\" name=\"f1\" spellcheck=\"false\" rows=\"15\">" << str1 << "</textarea>"
	    << "</td><td>&nbsp;</td><td>"
 	    << "<label for=\"f2\">Annotations</label><br>"
	    << "<textarea style=\"width: 100%; max-width: 100%; resize: none; font-family: Courier New\" id=\"f2\" name=\"f2\" spellcheck=\"false\" rows=\"15\">" << str2 << "</textarea>"
    	    << "</td></tr>"

	    << "<tr><td>"    
    	    << "<label for=\"f3\">Optional channel aliases <em>(from : to)</em></label><br>"
	    << "<textarea style=\"width: 100%; max-width: 100%; resize: none; font-family: Courier New\" id=\"f3\" name=\"f3\" spellcheck=\"false\" rows=\"10\">"<<str3<<"</textarea>"
	    << "</td><td>&nbsp;</td><td>"
 	    << "<label for=\"f4\">Optional annotation aliases <em>(from : to)</em></label><br>"
	    << "<textarea style=\"width: 100%; max-width: 100%; resize: none; font-family: Courier New\" id=\"f4\" name=\"f4\" spellcheck=\"false\" rows=\"10\">"<<str4<<"</textarea>"
    	    << "</td></tr>"

	    << "</table>"
	    << "<input type=\"submit\" value=\"Submit\">"
    	    << "<input type=\"reset\" value=\"Reset\">"
	    << "</form>"; 

}

void output_page( const std::map<std::string,std::string> & vars )
{

  std::map<std::string,std::string>::const_iterator ii = vars.begin();

  if ( vars.find( "f1" ) == vars.end()
       || vars.find( "f2" ) == vars.end()
       || vars.find( "f3" ) == vars.end()
       || vars.find( "f4" ) == vars.end() )
    return;
  
  //
  // make call to cgi-mapper;
  //
  
  std::vector<std::string> tok;
    
  //
  // added NAP remappings
  //

  tok.push_back( "amap=harm.annots" );
  tok.push_back( "cs-harm=harm.canonical.sigs" );
  tok.push_back( "cs-base=base.canonical.sigs" );

  std::string str1, str2, str3, str4;
  
  // channels 'f1' ... 'allow | chars , i.e. from Luna output
  std::vector<std::string> ctok = Helper::quoted_parse( vars.find( "f1" )->second  , " \n\r" );
  for (int i=0; i<ctok.size(); i++)
    {
      std::string v = Helper::trim( ctok[i] );
      if ( v != "|" )
	{
	  tok.push_back( "c=" + v );
	  if ( str1 != "" ) str1 += "|";
	  str1 += v;
	}
      
    }

  
  // annots 'f2' : also ignore '|' and '(...'  to make parsing Luna output easier...
  std::vector<std::string> atok = Helper::quoted_parse( vars.find( "f2" )->second , " \n\r" );
  for (int i=0; i<atok.size(); i++)
    {
      std::string v = Helper::trim( atok[i] );
      if ( v != "|" && v.substr(0,1) != "(" )
	{
	  tok.push_back( "a=" + v );
	  if ( str2 != "" ) str2 += "|";
          str2 += v;
	}
      
    }
  

  // user-specified channel mappings 'f3' (row delimited) 
  ctok = Helper::quoted_parse( vars.find( "f3" )->second , "\n\r" );
  for (int i=0; i<ctok.size(); i++)
    {
      std::vector<std::string> ttok = Helper::quoted_parse( ctok[i] , " " );
      if ( ttok.size() == 2 )
	{	  
	  tok.push_back( "alias=\"" + ttok[1] + "\"|\"" + ttok[0] + "\"" );
	  if ( str3 != "" ) str3 += "|";
	  str3 +=  ttok[0] + "=" + ttok[1];
	}

    }

  
  // user-specified annotatopn mappings 'f4' (row delimited) 
  atok = Helper::quoted_parse( vars.find( "f4" )->second , "\n\r" );
  for (int i=0; i<atok.size(); i++)
    {
      std::vector<std::string> ttok = Helper::quoted_parse( atok[i] , " " );
      if ( ttok.size() == 2 )
	{
	  tok.push_back( "remap=\"" + ttok[1] + "\"|\"" + ttok[0] + "\"" );
	  if ( str4 != "" ) str4 += "|";
	  str4 +=  ttok[0] + "=" + ttok[1];
	}
    }

  //
  // link back to main page
  //

  std::cout << "<table width=100%><tr><td style=\"text-align:right\">";
  std::cout << "<a href=\"/cgi-bin/cgi-mapper?"
	    << "f1=" << str1 << "&"
	    << "f2=" << str2 << "&"
	    << "f3=" << str3 << "&"
	    << "f4=" << str4 << "&"
	    << "inp=1" ;
  std::cout << "\">(return)</a>"
	    << "</td></tr></table>";
  
  
  //
  // call Luna mapper function
  //

  const bool html_mode = true;

  Helper::channel_annot_mapper( tok , html_mode );
  
  
  //
  // add another link to get back to front page
  //

  
  std::cout << "<hr><table width=100%><tr><td style=\"text-align:right\">";

  std::cout << "<a href=\"/cgi-bin/cgi-mapper?"
	    << "f1=" << str1 << "&"
	    << "f2=" << str2 << "&"
	    << "f3=" << str3 << "&"
	    << "f4=" << str4 << "&"
	    << "inp=1" ;
  std::cout << "\">(return)</a>"
	    << "</td></tr></table>";

    
}

