
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

#include <iostream>
#include <sstream>
#include <ctime>
#include <string>
#include <iomanip>

#include "defs/defs.h"

class logger_t
{

 private:

  const std::string _log_header;

  std::ostream& _out_stream;
  
  std::stringstream ss;
  
  bool         is_off;
  
 public:
  
 logger_t( const std::string & log_header  ,
	  std::ostream& out_stream = std::cerr)
   : _log_header( log_header ) , _out_stream( out_stream ) 
  {
    is_off = false;      
  }
  
  void flush() { _out_stream.flush(); } 

  void off() { flush(); is_off = true; } 

  void banner( const std::string & v , const std::string & bd ) 
  {

    if ( is_off || globals::silent ) return;
    
    // initialize log with this message
    
    time_t rawtime;
    time (&rawtime);
    struct tm * timeinfo = localtime (&rawtime);
    
    char BUFFER[50];
    strftime(BUFFER, sizeof(BUFFER), "%d-%b-%Y %H:%M:%S", timeinfo); 

    _out_stream << "===================================================================" << "\n"
		<< _log_header
		<< " | " << v << ", " << bd << " | starting " << BUFFER  << " +++\n"
		<< "===================================================================" << std::endl;
  }

                         
   
  ~logger_t()
    {
      
      if ( is_off || globals::silent ) return;
      
      if ( ! globals::silent ) 
	{
	  
	  time_t rawtime;
	  time (&rawtime);
	  struct tm * timeinfo = localtime (&rawtime);
	  
	  char BUFFER[50];
	  strftime(BUFFER, sizeof(BUFFER), "%d-%b-%Y %H:%M:%S", timeinfo); 
	  
	  _out_stream << "-------------------------------------------------------------------"
		      << "\n"
		      << "+++ luna | finishing "
		      << BUFFER
		      << "                       +++\n"
		      << "==================================================================="
		      << std::endl;
	}

    }


  void warning( const std::string & msg )
  {
    if ( is_off ) return ;

    if ( globals::logger_function )
      (*globals::logger_function)( " ** warning: " + msg + " **" );
    else if ( globals::Rmode && globals::Rdisp )
      ss << " ** warning: " << msg << " ** " << std::endl;
    else
      _out_stream << " ** warning: " << msg << " ** " << std::endl;
  }
  
  
  template<typename T>           
    logger_t& operator<< (const T& data) 
    {
      if ( is_off ) return *this;      
      
      if ( globals::logger_function )
	{
	  std::stringstream ss1;
	  ss1 << data;
	  (*globals::logger_function)( ss1.str() );
	}
      else if ( globals::Rmode && globals::Rdisp )
	ss << data;
      else if ( ! globals::silent ) 
	_out_stream << data;
      
      return *this;

    }
  

  std::string print_buffer() 
    {      
      std::string retval = ss.str();
      ss.str(std::string());
      return retval;
    }
  

};


#endif
