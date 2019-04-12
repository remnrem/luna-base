
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

#include "edfz/edfz.h"

#include <iostream>
#include "bgzf.h"
#include <cstdlib>
#include <vector>
#include <map>
#include "helper/helper.h"


// int test()
// {

//   edfz_t edfz;

  
//   bool okay = edfz.open_for_writing( "my.edfz" );

//   std::map<int64_t,int> offsets;

//   // read files from STDIN
  
//   while ( !std::cin.eof() ) 
//     {
//       std::string line;
      
//       std::getline( std::cin , line );

//       if ( std::cin.eof() ) break;

//       // add back
//       line += "\n";

//       int64_t offset = edfz.tell();
      
//       offsets[ offset ] = line.size();

//       edfz.write( (byte_t*)line.c_str() , line.size() );
	
      
//     }
  
//   edfz.close();
    
  
//   //
//   // Now, track through backwards
//   //

//   edfz_t reader;
  
//   reader.open_for_reading( "my.edfz" );

//   std::cout << "now going to read\n";

//   std::map<int64_t,int>::const_iterator ii = offsets.begin();
//   while ( ii != offsets.end() )
//     {

//       int64_t offset = ii->first;

//       int n = ii->second;
      
//       std::cout << offset << " " << n << "\n";

//       std::vector<byte_t> buffer( n );
      
//       byte_t * p = &(buffer[0]);
      
//       reader.read( p , n );

//       std::cout << "[";
//       for (int i=0;i<n;i++) std::cout << (char)buffer[i];
//       std::cout << "]\n\n";

//       ii++;
//     }

//   reader.close();
    

  
//   return 0;

// }
