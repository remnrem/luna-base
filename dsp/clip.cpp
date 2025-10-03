
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

#include "dsp/clip.h"

#include "edf/edf.h"
#include "edf/slice.h"

#include "param.h"

#include "helper/helper.h"
#include "helper/logger.h"

extern logger_t logger;

void dsptools::clip( edf_t & edf , param_t & param )
{

  //
  // get signals
  //
  
  std::string signal_label = param.requires( "sig" );
  const bool no_annotations = true;
  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );
  const int ns = signals.size();

  //
  // clipping parameters
  //

  bool abs_min = param.has( "lwr" );
  bool abs_max = param.has( "upr" );
  bool pct_min = param.has( "lwr-pct" );
  bool pct_max = param.has( "upr-pct" );

  double abs_min_th = abs_min ? param.requires_dbl( "lwr" ) : 0 ;
  double abs_max_th = abs_max ? param.requires_dbl( "upr" ) : 0 ;
  double pct_min_th = pct_min ? param.requires_dbl( "lwr-pct" ) : 0 ;
  double pct_max_th = pct_max ? param.requires_dbl( "upr-pct" ) : 0 ;

  if ( pct_min && ( pct_min_th <= 0 || pct_min_th >= 1 ) )
    Helper::halt( "pct_min_th must be between 0 and 1" ); 

  if ( pct_max && ( pct_max_th <= 0 || pct_max_th >= 1 ) )
    Helper::halt( "pct_max_th must be between 0 and 1" ); 
  
  //
  // nothing to do 
  //

  if ( ! ( abs_min_th || abs_max_th || pct_min_th || pct_max_th ) )
    return;
  
  //
  // process data 
  //

  logger << "  clipping signals:";
  
  for (int s=0; s<ns; s++)
    {

      logger << " " << signals.label(s) ;

      slice_t slice( edf , signals(s) , edf.timeline.wholetrace() );
      
      std::vector<double> * d = slice.nonconst_pdata();
      
      const int n = d->size();
      
      // absolute min threshold
      if ( abs_min ) 
	for (int i=0; i<n; i++)
	  if (  (*d)[i] < abs_min_th ) (*d)[i] = abs_min_th;
      
      // absolute max threshold
      if ( abs_max ) 
	for (int i=0; i<n; i++)
	  if (  (*d)[i] > abs_max_th ) (*d)[i] = abs_max_th;
      
      
      // pct min threshold
      if ( pct_min )
	{
	  const double p = MiscMath::percentile( *d , pct_min_th );	  
	  for (int i=0; i<n; i++)
	    if (  (*d)[i] < p ) (*d)[i] = p;
	}
      
      // pct max threshold
      if ( pct_max )
	{
	  const double p = MiscMath::percentile( *d , pct_max_th );	      
	  for (int i=0; i<n; i++)
	    if (  (*d)[i] > p ) (*d)[i] = p;
	}
      
      edf.update_signal( signals(s) , d );
      
    }
  
  logger << "\n";
  
}


