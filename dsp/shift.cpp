
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

#include "dsp/shift.h"

#include "param.h"
#include "edf/edf.h"

#include <iostream>
#include <cmath>

void dsptools::shift( edf_t & edf , param_t & param ) 
{

  std::string signal_label = param.value( "sig" );

  signal_list_t signals = edf.header.signal_list( signal_label );  

  const int ns = signals.size();
  
  const int nt = param.requires_int( "sp" );
  
  bool nowrap = param.has( "no-wrap" );

  for ( int s = 0 ; s < ns ; s++ )
    {
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;	  

      logger << "  shifting " << signals.label(s) << " by " << nt << " sample points";
      if ( nowrap ) logger << " (no wrapping)\n";
      else logger << " (wrapping)\n";
      
      edf.shift( signals(s) , nt , ! nowrap );

    }
  
}

