
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

#include "sstore.h"
#include <iostream>
#include "helper/helper.h"

//
// loadints, a really simple interval loader for a sstore_t 
// 

int main(int argc , char ** argv )
{

  if ( argc != 3 ) 
    { 
      std::cerr << "usage: ./loadints {filename} {label} < input\n"
		<< " where input is the interval ranges (seconds)"
		<< "\n";          
	std::exit(1); 
    } 
  

  // 
  // Input format, tab-delimited
  //
  
  // interval :  START STOP

  std::string filename = argv[1];
  std::string label    = argv[2];
  
  //
  // Open/create sstore_t
  //
  
  sstore_t ss( filename );
  
  ss.drop_index();
  
  while ( ! std::cin.eof() ) 
    {
      double a,b,;
      std::cin >> a >> b;
      if ( std::cin.eof() ) break;
      
      // only label (as name) : no value, or channel/level stratifiers
      ss.insert_interval( a , b , label , "." , NULL , NULL );
      
      // next row of input
    }
  
  ss.index();

  std::exit(0);
}
