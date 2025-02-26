
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


#ifndef __PSGTOOLS_MAIN_H__
#define __PSGTOOLS_MAIN_H__

#include <string>
#include <map>
#include <set>
#include <sstream>
#include <iostream>
#include <new>
#include <fstream>

#include "miscmath/crandom.h"

class param_t;
class cmd_t;

// misc helper: command logger 
std::string log_commands( int argc , char ** argv );

// misc helper: sample-list slicer
bool luna_helper_sl_slicer( const std::string & f , int n , int m , int * s1 , int * s2 );

// misc helper: evaluate comannd-line methods
cmdline_proc_t parse_cmdline( int argc , char ** argv , int * );

// misc helper: execture special cmd-line ops
void exec_cmdline_procs( cmdline_proc_t & cmdline , int argc , char ** argv, int param_from_command_line );

// misc helper: read and parse an @include file
void include_param_file( const std::string & paramfile );

// misc helper: manage memory resource issues
void NoMem();

void proc_eval_tester( const bool );

void process_edfs(cmd_t&);

// misc helper: build global params from cmdline
void build_param( param_t * , int argc , char** argv , int );

// misc helper: build global params from stdin
void build_param_from_stdin( param_t * );

// misc helper: return Luna version
std::string luna_base_version();

#endif
