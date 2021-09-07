
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

#include "utils/cgi-utils.h"

#include <iostream>


int main( int argc , char ** argv )
{
  
  html_write_headers( "Test One" );

  std::map<std::string,std::string> vars = fetch_cgi();

  std::map<std::string,std::string>::const_iterator ii = vars.begin();
  while ( ii != vars.end() )
    {
      std::cout << ii->first << " = " << ii->second << "</p>";
      ++ii;
    }
  

  html_write_footer();

}
  
