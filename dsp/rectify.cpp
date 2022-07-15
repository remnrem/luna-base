
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

#include "dsp/rectify.h"

#include "edf/edf.h"
#include "edf/slice.h"

#include "helper/helper.h"
#include "helper/logger.h"

extern logger_t logger;

void dsptools::rectify( edf_t & edf , param_t & param )
{

  //
  // get signals
  //
  
  std::string signal_label = param.requires( "sig" );
  const bool no_annotations = true;
  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );
  const int ns = signals.size();

  //
  // process data 
  //

  logger << "  rectifying signals:";
  
  for (int s=0; s<ns; s++)
    {

      logger << " " << signals.label(s) ;
      
      slice_t slice( edf , signals(s) , edf.timeline.wholetrace() );

      std::vector<double> * d = slice.nonconst_pdata();

      const int n = d->size();
      
      for (int i=0; i<n; i++)
	(*d)[i] = fabs( (*d)[i] ) ;
      
      edf.update_signal( signals(s) , d );
      
    }

  logger << "\n";
  
}


