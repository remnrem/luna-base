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


#include "spindles/spindles.h"

#include "helper/logger.h"

#include "db/db.h"

extern logger_t logger;

extern writer_t writer;

void spindle_so_coupling( edf_t & edf , param_t & param )
{

  const std::string spindle_cache = param.requires( "spindles" );
  const std::string so_cache = param.requires( "so" );

  //
  // Options
  //

  // spindle/SWA or spindle/SO coupling?
  bool all_spindles = param.has( "all-spindles" );

  // use randomization (shuffled) surrogate time-series?
  int nreps = param.has( "nreps" ) ? param.requires_int( "nreps" ) : 0;

  // shuffle only within epoch?
  bool stratify_by_so_phase_bin = param.has ( "stratify-by-phase" );

  // Within-epoch permutation (default, unless this given)
  bool eperm = ! param.has( "perm-whole-trace" );
  

  //
  // report what we're doing here:
  //

  logger << "  spindle/SO coupling\n";
    
  
  //
  // retrieve pre-calculated values from the cache
  //
  
  
  //
  // determine unique channels
  //

  
  //
  // iterate over each class of SP and SO, but only within channel
  //

  
  
}
