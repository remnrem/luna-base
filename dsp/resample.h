
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

#ifndef __RESAMPLE_H__
#define __RESAMPLE_H__

#include <vector>
#include <string>

#include "libsamplerate/samplerate.h"


// interface to SRC libsamplerate (which must be installed on system)
// http://www.mega-nerd.com/SRC/

struct edf_t;
struct param_t;

namespace dsptools 
{

  void resample_channel( edf_t & , param_t & );

  void resample_channel_zoh( edf_t & , param_t & );

  void resample_channel( edf_t & , const int , const double , const int converter = SRC_SINC_FASTEST );

  std::vector<double> resample( const std::vector<double> * d , double sr1 , double sr2 , int converter = SRC_SINC_FASTEST );

  int converter( const std::string & m );

  std::string converter( int m );

  void fix_sampling( edf_t & , param_t & );
  
}

  
#endif
