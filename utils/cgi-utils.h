
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

#ifndef __CGI_UTILS_H__
#define __CGI_UTILS_H__


/***************** start of the getcgivars() module **********************/

/*************************************************************************/
/**                                                                     **/
/**     getcgivars.c-- routine to read CGI input variables into an      **/
/**         array of strings.                                           **/
/**                                                                     **/
/**     Written in 1996 by James Marshall, james@jmarshall.com, except  **/
/**     that the x2c() and unescape_url() routines were lifted directly **/
/**     from NCSA's sample program util.c, packaged with their HTTPD.   **/
/**                                                                     **/
/**     For the latest, see http://www.jmarshall.com/easy/cgi/ .        **/
/**                                                                     **/
/*************************************************************************/

#include <cstdio>
#include <cstring>
#include <cstdlib>

/** Convert a two-char hex string into the char it represents. **/
char x2c(char *what) {
  char digit;  
  digit = (what[0] >= 'A' ? ((what[0] & 0xdf) - 'A')+10 : (what[0] - '0'));
  digit *= 16;
  digit += (what[1] >= 'A' ? ((what[1] & 0xdf) - 'A')+10 : (what[1] - '0'));
  return(digit);
}

/** Reduce any %xx escape sequences to the characters they represent. **/
void unescape_url(char *url) {
  int i,j;
  
  for(i=0,j=0; url[j]; ++i,++j) {
    if((url[i] = url[j]) == '%') {
      url[i] = x2c(&url[j+1]) ;
      j+= 2 ;
    }
  }
  url[i] = '\0' ;
}


/** Read the CGI input and place all name/val pairs into list.        **/
/** Returns list containing name1, value1, name2, value2, ... , NULL  **/

char **getcgivars() {
    int i ;
    char *request_method ;
    int content_length;
    char *cgiinput ;
    char **cgivars ;
    char **pairlist ;
    int paircount ;
    char *nvpair ;
    char *eqpos ;

    /** Depending on the request method, read all CGI input into cgiinput. **/
    request_method= getenv("REQUEST_METHOD") ;
    
    if (!strcmp(request_method, "GET") || !strcmp(request_method, "HEAD") ) {
      /* Some servers apparently don't provide QUERY_STRING if it's empty, */
      /*   so avoid strdup()'ing a NULL pointer here.                      */
      char *qs ;
      qs= getenv("QUERY_STRING") ;
      cgiinput= strdup(qs  ? qs  : "") ;
    }
    else if (!strcmp(request_method, "POST")) {
      /* strcasecmp() is not supported in Windows-- use strcmpi() instead */
      if ( strcasecmp(getenv("CONTENT_TYPE"), "application/x-www-form-urlencoded")) {
	printf("Content-Type: text/plain\n\n") ;
	printf("getcgivars(): Unsupported Content-Type.\n") ;
	exit(1) ;
      }
      if ( !(content_length = atoi(getenv("CONTENT_LENGTH"))) ) {
	printf("Content-Type: text/plain\n\n") ;
	printf("getcgivars(): No Content-Length was sent with the POST request.\n") ;
	exit(1) ;
      }
      if ( !(cgiinput= (char *) malloc(content_length+1)) ) {
	printf("Content-Type: text/plain\n\n") ;
	printf("getcgivars(): Couldn't malloc for cgiinput.\n") ;
	exit(1) ;
      }
      if (!fread(cgiinput, content_length, 1, stdin)) {
	printf("Content-Type: text/plain\n\n") ;
	printf("getcgivars(): Couldn't read CGI input from STDIN.\n") ;
	exit(1) ;
      }
      cgiinput[content_length]='\0' ;
    }
    else {
      printf("Content-Type: text/plain\n\n") ;
      printf("getcgivars(): Unsupported REQUEST_METHOD.\n") ;
      exit(1) ;
    }
    
    /** Change all plusses back to spaces. **/
    for (i=0; cgiinput[i]; i++)
      if (cgiinput[i] == '+') cgiinput[i] = ' ' ;
    
    /** First, split on "&" and ";" to extract the name-value pairs into **/
    /**   pairlist.                                                      **/
    pairlist= (char **) malloc(256*sizeof(char **)) ;
    paircount= 0 ;
    nvpair= strtok(cgiinput, "&;") ;
    while (nvpair) {
      pairlist[paircount++]= strdup(nvpair) ;
      if (!(paircount%256))
	pairlist= (char **) realloc(pairlist,(paircount+256)*sizeof(char **)) ;
      nvpair= strtok(NULL, "&;") ;
    }
    pairlist[paircount]= 0 ;    /* terminate the list with NULL */
    
    /** Then, from the list of pairs, extract the names and values. **/
    cgivars= (char **) malloc((paircount*2+1)*sizeof(char **)) ;
    for (i= 0; i<paircount; i++) {
      if ( ( eqpos = strchr(pairlist[i], '=') ) ) {
	*eqpos= '\0' ;
	unescape_url(cgivars[i*2+1]= strdup(eqpos+1)) ;
      } else {
	unescape_url(cgivars[i*2+1]= strdup("")) ;
      }
      unescape_url(cgivars[i*2]= strdup(pairlist[i])) ;
    }
    cgivars[paircount*2]= 0 ;   /* terminate the list with NULL */
    
    /** Free anything that needs to be freed. **/
    free(cgiinput) ;
    for (i=0; pairlist[i]; i++) free(pairlist[i]) ;
    free(pairlist) ;
    
    /** Return the list of name-value strings. **/
    return cgivars ;
    
}

/***************** end of the getcgivars() module ********************/



#include <map>
#include <iostream>
#include <string>


void html_write_headers( const std::string & title )
{
  std::cout << "Content-type: text/html\n\n"
    "<!DOCTYPE html>"
    "<html lang=\"en\">"
    "<head><meta charset=\"utf-8\"><title>" << title << "</title>"
    "<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\" />";
  
  std::cout << "<html><body>";
  
}

void html_write_footer( )
{
  std::cout << "<br><hr></body></html>";
}


std::map<std::string,std::string> fetch_cgi() {
  std::map<std::string,std::string> vals;  
  char ** cgivars = getcgivars();  
  for (int i=0; cgivars[i]; i+= 2)
    vals[ cgivars[i] ] = cgivars[i+1];  
  // free CGI variables  
  for (int i=0; cgivars[i]; i++) free(cgivars[i]) ;
  free(cgivars) ;  
  return vals;
}

#include <iostream>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>

std::string exec_system( const std::string & cmd ) {
  
  std::array<char, 128> buffer;
  
  std::string result;
  
  std::unique_ptr<FILE, decltype(&pclose)> pipe(popen( cmd.c_str() , "r"), pclose);
  
  if ( ! pipe ) {
    throw std::runtime_error("popen() failed!");
  }
  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
    result += buffer.data();
  }
  return result;
}

#endif
