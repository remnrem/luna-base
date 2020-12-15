
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

#include "timeline/cache.h"

#include "db/db.h"

extern writer_t writer;

void ctest()
{
  
  writer.level( "L1" , "F1" );
  
  writer.level( 123 , "FFE" );
  writer.epoch( 222 );
  
  cache_t<double> cache( "my1" );

  ckey_t ckey1( "y" , writer.faclvl() ) ;
  ckey_t ckey2( "z" , writer.faclvl() ) ;

  std::vector<double> y( 10 , 22 );
  std::vector<double> z( 10 , 23 );
  
  cache.add( ckey1 , y );
  cache.add( ckey2 , z );
  
  writer.unlevel();

  std::cout << cache.print();
  
}
