
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

// log utility initially based on: https://github.com/Manu343726/Cpp11CustomLogClass

#ifndef __LOGGER_H__
#define	__LOGGER_H__

#include <type_traits>
#include <iostream>
#include <chrono>
#include <ctime>
#include <string>
#include <iomanip>

#include "defs/defs.h"

using std::chrono::system_clock;

class logger_t
{

 private:

  std::ostream& _out_stream;
  
  bool         _next_is_begin;
  
  const std::string _log_header;

  using endl_type = std::ostream&(std::ostream&);

  //This is the key: std::endl is a template function, and this is the
  //signature of that function (For std::ostream).
  
 public:
  
  // Constructor: User passes a custom log header and output stream, or uses defaults.

   logger_t(const std::string& log_header  ,
	  std::ostream& out_stream = std::cerr)
   : _log_header( log_header ) , _out_stream( out_stream ) , _next_is_begin( true )
    {

      // initialize log with this message
      
      auto now        = std::chrono::system_clock::now();
      auto now_time_t = std::chrono::system_clock::to_time_t( now ); 
      auto now_tm     = std::localtime( &now_time_t ); 

      if ( ! globals::silent ) 
	_out_stream << _log_header
		    << " | starting process "
		    << std::put_time( now_tm, "%Y-%m-%d %H:%M:%S")		  
		    << std::endl;
      
    }



  ~logger_t()
    {
      
      // initialize log with this message
      
      auto now        = std::chrono::system_clock::now();
      auto now_time_t = std::chrono::system_clock::to_time_t( now ); 
      auto now_tm     = std::localtime( &now_time_t ); 

      if ( ! globals::silent ) 
	_out_stream << "-------------------------------------------------------------------"
		    << std::endl 
		    << "+++ luna | finishing process "
		    << std::put_time( now_tm, "%Y-%m-%d %H:%M:%S")
		    << std::endl
		    << "==================================================================="
		    << std::endl;
      
    }

  
  
  // Overload for std::endl only:
  logger_t& operator<<(endl_type endl)
    {       
      _next_is_begin = true; 

      if ( ! globals::silent ) 
	_out_stream << endl; 

      return *this; 
    }
  
  //Overload for anything else:

  template<typename T>           
    logger_t& operator<< (const T& data) 
    {

      if ( ! globals::silent ) 
	_out_stream << data;
      

      /* auto now        = std::chrono::system_clock::now(); */
      /* auto now_time_t = std::chrono::system_clock::to_time_t( now ); //Uhhg, C APIs... */
      /* auto now_tm     = std::localtime( &now_time_t ); //More uhhg, C style...  */

      /* auto t = std::time(nullptr); */
      /* auto tm = *std::localtime(&t); */

      /* << "(" */
      /* << ( now_tm->tm_hour < 10 ? "0" : "" ) << now_tm->tm_hour << ":" */
      /* << ( now_tm->tm_min < 10 ? "0" : "" ) << now_tm->tm_min << ":" */
      /* << ( now_tm->tm_sec < 10 ? "0" : "" ) << now_tm->tm_sec << "): " */

      
      /* if( _next_is_begin ) */
      /* 	_out_stream << _log_header */
      /* 		    << " (" */
      /* 		    << std::put_time( now_tm, "%Y-%m-%d %H:%M:%S") */
      /* 		    << "): " */
      /* 		    << data; */
      /* else */
      /* 	_out_stream << data; */
      
      /* _next_is_begin = false; */

      return *this;

    }
};


#endif
