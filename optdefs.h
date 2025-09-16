
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

#ifndef __LUNA_OPTDEFS_H__
#define __LUNA_OPTDEFS_H__

#include <string>
#include <map>
#include "eval.h"
#include "helper/helper.h"
#include <vector>
#include <set>
#include <sstream>

//
// optdefs_t handles a) knowledge about Luna spectial variables
//

enum opt_type_t {
  OPT_FLAG_T ,
  OPT_BOOL_T ,
  OPT_INT_T ,
  OPT_NUM_T ,
  OPT_NUM_INTERVAL_T , 
  OPT_STR_T , 
  OPT_FILE_T ,
  OPT_PATH_T ,
  OPT_INTVEC_T ,
  OPT_NUMVEC_T ,
  OPT_STRVEC_T ,
  OPT_CHAR_T , 
  OPT_TIME_T ,
  OPT_DATE_T , 
  OPT_SPECIAL_T ,
  OPT_UNDEFINED_T 
};


class optdefs_t 
{
  
 public:
  
  optdefs_t();
  
  void init();
  
  void add( const std::string & domain, const std::string & opt , const opt_type_t & type , const std::string & desc );

  bool has( const std::string & opt ) const
  {
    return odesc.find( opt ) == odesc.end(); 
  }
  
  opt_type_t get_type( const std::string & d ) const;
  
  std::string get_desc( const std::string & d ) const;

  std::vector<std::string> get_domains( const std::string & ) const
  {
    return domains;
  }

  std::vector<std::string> get_opts( const std::string & ) const; 

  static std::string type( opt_type_t t )
  {
    if ( t == OPT_FLAG_T ) return "flag";
    else if ( t == OPT_BOOL_T ) return "true/false";
    else if ( t == OPT_INT_T ) return "integer";
    else if ( t == OPT_NUM_T ) return "numeric";
    else if ( t == OPT_STR_T ) return "text";
    else if ( t == OPT_FILE_T ) return "file";
    else if ( t == OPT_PATH_T ) return "path";
    else if ( t == OPT_INTVEC_T ) return "integer-list";
    else if ( t == OPT_NUMVEC_T ) return "numeric-list";
    else if ( t == OPT_STRVEC_T ) return "text-list";
    else if ( t == OPT_SPECIAL_T ) return "special";
    else if ( t == OPT_CHAR_T ) return "char";
    else if ( t == OPT_TIME_T ) return "time";
    else if ( t == OPT_DATE_T ) return "date";
    else return "undefined";
  }

private:
  
  std::vector<std::string> domains;

  std::map<std::string,std::vector<std::string> > domain2opt;
  
  std::map<std::string,std::string> odesc;

  std::map<std::string,opt_type_t> otype;
  
};


#endif


