
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

#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <cstring>
#include <fstream>


struct reg_t {  

  reg_t( double a, double b) : start(a), stop(b) { }
  
  double start;
  double stop;

  bool operator<( const reg_t & rhs ) const 
  {
    if ( start < rhs.start ) return true;
    if ( start > rhs.start ) return false;
    return stop < rhs.stop;
  }
};


int main(int argc , char ** argv )
{
  
  // regional
  //     - list of regions
  //     - 0 or 1 for includ if in, versus exclude if in
  //     - text file
  //     - 1 or 2
  //     - if 1 :   VAR and then +/2 window (numeric)
  //     - if 2 :   VAR1 and VAR2
  
  // output:
  //  lines from text file that match  
  
  if ( argc != 7 )
    {
      std::cerr << " error: usage \n"
		<< "   regional <regions> <0|1> <input> <1|2> <VAR|VAR1> <W|VAR2> \n";
      std::exit(1);
    }
  
  const std::string region_file = argv[1] ;

  const bool include =  strcmp( argv[2] , "1" ) == 0 ;
  
  const std::string input_file = argv[3];
  
  const bool single_point = strcmp( argv[4] , "1" ) == 0 ;

  const std::string var1 = argv[5];

  const std::string var2 = single_point ? "" : argv[6];

  const double w = single_point ? atof( argv[6] ) : 0 ;

  std::cerr << " args\n"
	    << "  regions = " << region_file << "\n"
	    << "  include = " << ( include ? "yes" : "no" ) << "\n"
	    << "  input   = " << input_file << "\n"
	    << "  regions = " << ( single_point ? "single-point" : "start-stop" ) << "\n";

  if ( single_point )
    std::cout << "  center var  = " << var1 << "\n"
	      << "  window  = " << w << "\n";
  else
    std::cout << "  start var  = " << var1 << "\n"
	      << "  stop var  = " << var2 << "\n";


  // read regions
  
  std::exit(0);
}
