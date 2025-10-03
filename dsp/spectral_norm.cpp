
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

#include "spectral_norm.h"

#include <string>

#include "edf/edf.h"
#include "edf/slice.h"

#include "param.h"

#include "helper/logger.h"

extern logger_t logger;


// implements https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2235870/

void dsptools::norm_1overf( edf_t & edf , param_t & param )
{
  
  std::string signal_label = param.requires( "sig" );
  signal_list_t signals = edf.header.signal_list( signal_label );    
  const int ns = signals.size();

  for (int s=0;s<ns;s++)
    {
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
      
      const double fs = edf.header.sampling_freq( signals(s) ) ; 
      
      logger << "  1/f normalizing " <<  signals.label(s)  << "(Fs=" << fs << ")\n";

      // get data 

      interval_t interval = edf.timeline.wholetrace();
      
      slice_t slice( edf , signals(s) , interval );    

      const std::vector<double> * sig  = slice.pdata();

      std::vector<double> nsig = dsptools::norm_1f( *sig , fs );
      
      edf.update_signal( signals(s) , &nsig );      
      
    }

}

std::vector<double> dsptools::norm_1f( const std::vector<double> & x , double fs )
{
  const int n = x.size();
  const double dt = 1.0 / (double) fs;
  std::vector<double> y( n , 0 );
  for (int i=1;i<n;i++) y[i] = ( x[i] - x[i-1] ) / dt;
  return y;
}

