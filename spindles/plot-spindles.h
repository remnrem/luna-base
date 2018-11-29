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


#ifndef __SPINDLE_PLOTS_H__
#define __SPINDLE_PLOTS_H__

void draw_spindles( edf_t & edf , 
		    param_t & param , 
		    const std::string & filename ,
		    int s , 
		    const std::vector<spindle_t> & spindles , 		    
		    const std::map<uint64_t,double> * avgmap ) ;

void draw_mspindles( edf_t & edf , 
		     param_t & param , 
		     const std::string & filename ,
		     std::vector<int> s , 
		     const std::vector<mspindle_t> & spindles );




#endif
