
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

#ifndef __ZFILE_H__
#define __ZFILE_H__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <set>


#ifdef __WIN
#include <direct>
#else
#include <sys/stat.h>
#endif

#include <sys/stat.h>
#include <sys/types.h>
#include <cerrno>

#include <cstdio>
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#include "zfstream.h"

#include <wordexp.h>


// wrappers for input and output streams

class InFile : public gzifstream {
 public:

  InFile() { } 
  
 InFile( const std::string & n , 
	 ios_base::openmode mode = ios_base::in ) 
   : gzifstream(n.c_str(), mode) { Helper::checkFileExists(n); } 
  
  std::string readLine();
  std::vector< std::string > tokenizeLine( const std::string & delim = PLINKSeq::DELIM() );
};

class OutFile : public gzofstream {
 public:
 OutFile( const std::string & n , 
	  ios_base::openmode mode = ios_base::out ) 
   : gzofstream(n.c_str(), mode) { }   
};

#endif
