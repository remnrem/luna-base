
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

#include "defs/defs.h"


class logger_t
{

 private:

  const std::string _log_header;

  std::ostream& _out_stream;
  
  bool         _next_is_begin;
 
  bool         is_off;

  using endl_type = std::ostream&(std::ostream&);


  //This is the key: std::endl is a template function, and this is the
  //signature of that function (For std::ostream).
  
 public:
  
  // Constructor: User passes a custom log header and output stream, or uses defaults.

   logger_t(const std::string& log_header  ,
	  std::ostream& out_stream = std::cerr)
   : _log_header( log_header ) , _out_stream( out_stream ) , _next_is_begin( true )
    {
      is_off = false;      
    }

  void flush() { _out_stream.flush(); } 

  void off() { flush(); is_off = true; } 

  void banner( const std::string & v , const std::string & bd ) 
  {

    if ( is_off || globals::silent ) return;

    // initialize log with this message
    
    auto now        = std::chrono::system_clock::now();
    auto now_time_t = std::chrono::system_clock::to_time_t( now ); 
    auto now_tm     = std::localtime( &now_time_t ); 
    
    _out_stream << "===================================================================" << std::endl
		<< _log_header
		<< " | " << v << ", " << bd 
		<< " | starting process "
		<< std::put_time( now_tm, "%Y-%m-%d %H:%M:%S")		  
		<< std::endl
		<< "===================================================================" << std::endl;
  }
  
  
  
  ~logger_t()
    {
      
      if ( is_off || globals::silent ) return;
      
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

  void warning( const std::string & msg )
  {
    if ( is_off ) return ;
    _out_stream << " ** warning: " << msg << " ** " << std::endl;
  }
  

  // Overload for std::endl only:
  logger_t& operator<<(endl_type endl)
    {       
      if ( is_off ) return *this;
      _next_is_begin = true; 

      if ( ! globals::silent ) 
	_out_stream << endl; 

      return *this; 
    }
  
  //Overload for anything else:

  template<typename T>           
    logger_t& operator<< (const T& data) 
    {
      if ( is_off ) return *this;      

      if ( ! globals::silent ) 
	_out_stream << data;

      return *this;

    }
};


#endif
