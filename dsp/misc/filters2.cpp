
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


#include "filters2.h"


std::vector<double> filter_t::filter1( const std::vector<double> & x )
{  
  // forward
  int i = 0;
  std::vector<double> y( x.size() );
  std::vector<double>::const_iterator fi = x.begin();
  while ( fi != x.end() ) { y[i++] = f1->do_sample( *fi ); ++fi; }
  
  // reverse
  i = 0;
  std::vector<double>::reverse_iterator ri = y.rbegin();
  while ( ri != y.rend() ) { *ri = f1->do_sample( *ri ); ++ri; }
  
  return y;
}
