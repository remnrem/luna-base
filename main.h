
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

class param_t;
class cmd_t;

void proc_dummy( const std::string & );
void proc_eval_tester( const bool );
void process_edfs(cmd_t&);
void list_cmds();

void build_param_from_cmdline( param_t * );
std::string luna_base_version();

#endif
