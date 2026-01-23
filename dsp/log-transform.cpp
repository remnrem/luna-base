
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

#include "dsp/log-transform.h"

#include "param.h"
#include "edf/edf.h"
#include "edf/slice.h"

#include "helper/helper.h"
#include "helper/logger.h"

extern logger_t logger;

void dsptools::log_transform( edf_t & edf , param_t & param )
{

  //
  // get signals
  //
  
  std::string signal_label = param.requires( "sig" );
  const bool no_annotations = true;
  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );
  const int ns = signals.size();

  //
  // options
  //

  const double eps_th = param.has( "eps" ) ? param.requires_dbl( "eps" ) : 0.01 ;  

  if ( eps_th <= 0 || eps_th > 0.2 ) Helper::halt( "eps must be 0 < th <= 0.2" );

  logger << "  using eps = 0.1 * " << eps_th << " percentile\n";
  
  //
  // process data 
  //

  logger << "  log-transforming signals:";
  
  for (int s=0; s<ns; s++)
    {

      logger << " " << signals.label(s) ;
      
      slice_t slice( edf , signals(s) , edf.timeline.wholetrace() );

      std::vector<double> d = *slice.nonconst_pdata();

      const int n = d.size();

      // check for negative
      int neg_cnt = 0;
      double max_neg = 0;
      for (int i=0; i<n; i++)
	if ( d[i] < 0 )
	  {
	    if ( d[i] < max_neg ) max_neg = d[i];
	    ++neg_cnt;
	    d[i] = 0;
	  }

      if ( neg_cnt )
	logger << "\n"
	       << "  *** warning - " << neg_cnt 
	       << " (of " << n << ") negative values found, clamping to 0.0\n"
	       << "  ***         - minimum value = " << max_neg << "\n";

      // get eps
      double eps_val = 0.1 * MiscMath::percentile( d , eps_th );
      if (!(eps_val > 0.0) || !std::isfinite(eps_val)) eps_val = 1e-12;

      // do transform
      for (int i=0; i<n; i++)
	d[i] = log( d[i] + eps_val ) ;
      
      edf.update_signal( signals(s) , &d );
      
    }

  logger << "\n";
  
}


