
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
#include "miscmath/miscmath.h"

#include "helper/helper.h"
#include "helper/logger.h"

#include <cmath>

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
  // options:
  //  offset=<c>  fixed additive constant before log (e.g. offset=1 → log(x+1))
  //  eps=<p>     data-driven offset: 0.1 × p-th percentile (default p=0.01, i.e. 1st pctile)
  //              ignored if offset= is given
  //

  const bool use_fixed_offset = param.has( "offset" );
  double fixed_offset = 0;

  double eps_th = 0.01;

  if ( use_fixed_offset )
    {
      fixed_offset = param.requires_dbl( "offset" );
      if ( fixed_offset < 0 )
	Helper::halt( "LOG offset must be >= 0" );
      logger << "  LOG: using fixed offset " << fixed_offset << " -> log(x + " << fixed_offset << ")\n";
    }
  else
    {
      eps_th = param.has( "eps" ) ? param.requires_dbl( "eps" ) : 0.01;
      if ( eps_th <= 0 || eps_th > 0.2 )
	Helper::halt( "LOG eps must be between 0 and 0.2 (as a proportion, e.g. 0.01 = 1st percentile)" );
      logger << "  LOG: using data-driven offset = 0.1 x "
	     << Helper::dbl2str( eps_th * 100.0 , 4 ) << "th percentile\n";
    }

  //
  // process signals
  //

  logger << "  LOG: transforming:";

  for (int s = 0; s < ns; s++)
    {
      logger << " " << signals.label(s);

      slice_t slice( edf , signals(s) , edf.timeline.wholetrace() );

      std::vector<double> d = *slice.nonconst_pdata();

      const int n = d.size();

      // clamp negatives to 0
      int neg_cnt = 0;
      double max_neg = 0;
      for (int i = 0; i < n; i++)
	if ( d[i] < 0 )
	  {
	    if ( d[i] < max_neg ) max_neg = d[i];
	    ++neg_cnt;
	    d[i] = 0;
	  }

      if ( neg_cnt )
	logger << "\n  *** warning - " << neg_cnt << " (of " << n
	       << ") negative values clamped to 0; minimum was " << max_neg << "\n";

      // determine offset
      double offset;
      if ( use_fixed_offset )
	{
	  offset = fixed_offset;
	}
      else
	{
	  offset = 0.1 * MiscMath::percentile( d , eps_th );
	  if ( ! ( offset > 0.0 ) || ! std::isfinite( offset ) )
	    {
	      offset = 1e-12;
	      logger << "\n  *** note: eps percentile was <= 0 or non-finite; using fallback offset = 1e-12\n";
	    }
	}

      // apply log(x + offset)
      for (int i = 0; i < n; i++)
	d[i] = std::log( d[i] + offset );

      edf.update_signal( signals(s) , &d );
    }

  logger << "\n";

}
