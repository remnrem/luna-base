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

#ifndef __SIGDYN_H__
#define __SIGDYN_H__

struct edf_t;
struct param_t;

// SIGDYN : generic temporal-dynamics profiling of any existing signal
//   0) simple stage/cycle/stability-stratified descriptive statistics
//      (N/MEAN/SD/MIN/MAX/RANGE)                          -- new local code
//   1) whole-recording trend/decile/NREM-cycle summary  -- reuses qdynam_t
//   2) anchor/peri-event windowed averaging              -- native engine
//
// This is a thin driver that calls existing, unmodified qdynam_t machinery
// for mode 1, plus its own local mode-0/mode-2 code (mode 2's peri-event
// engine is self-contained here, not a wrapper around dsp/tlock.{h,cpp});
// see dynamics/sigdyn.cpp for details.

namespace dsptools
{
  void sigdyn( edf_t & edf , param_t & param );
}

#endif
