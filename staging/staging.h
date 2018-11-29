
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


#ifndef __STAGING_H__
#define __STAGING_H__


#include <string>
#include <map>
#include <vector>

struct edf_t;


struct zratio_t
{
  void calc( edf_t & edf , const std::string & signal_label );
  std::vector<double> zr2;  // 2-second values
  std::vector<double> zr30; // averaged over 30-seconds
};

struct staging_t
{
  zratio_t zratio;
};

#endif

