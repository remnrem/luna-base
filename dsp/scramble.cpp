
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

#include "dsp/scramble.h"
#include "edf/edf.h"
#include "edf/slice.h"

#include "miscmath/crandom.h"
#include <iostream>
#include <cmath>

void dsptools::scramble( edf_t & edf , param_t & param ) 
{
  
  std::string signal_label = param.value( "sig" );
  
  signal_list_t signals = edf.header.signal_list( signal_label );  
  
  const int ns = signals.size();
  
  for ( int s = 0 ; s < ns ; s++ )
    {
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;	  
      
      const int slot = signals(s);
      
      logger << "  scrambling " << signals.label(s) << " completely (sample-by-sample randomization)\n";
      
      // get data : note, this ignores EDF discontinuities      
      slice_t slice( edf , slot , edf.timeline.wholetrace() );
      
      const std::vector<double> * d = slice.pdata();      
      const int np = d->size();
      
      // Yates-Fisher shuffle of indices
      std::vector<int> idx( np , 0 );
      CRandom::random_draw( idx );
      
      // new signal
      std::vector<double> d2( np , 0 );
      for (int i=0;i<np;i++)
	d2[ idx[i] ] = (*d)[i];
      
      // push back
      edf.update_signal( slot , &d2 );
      
    }
  
}

