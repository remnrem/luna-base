
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

#ifdef HAS_LGBM


// analog of POPS but for within-stage, between-person classification
// focus on binary group problem (but extension to other models should be trivial, e.g. QTs)

// same basic structure as POPS
//  - generation of *epoch-level* metrics, per EDF, and save (e.g. binary intermediates)
//    + we can use exactly the same machinery as POPS does?
//    ? what about micro-arch transients, e.g.
//    ? is epoch-level the right level of analysis?
//    -> ability to include other metrics: coherence / pairwise measures
//    -> generic use of cache mechanism to include new metrics?
//  - methods to include indiv-level covariates into the data

//  - collate all data, compute any 'level-2 metrics'
//    -> but perhaps less interested in temporal smoothing?

#include "lgbm/assoc.h"

assoc_t::assoc_t() {
  a = 1;
}


  
#endif

