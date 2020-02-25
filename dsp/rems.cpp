

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

#include "rems.h"

#include "intervals/intervals.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"

#include <vector>
#include <iostream>
#include <sstream>

extern writer_t writer;
extern logger_t logger;

void dsptools::rems( edf_t & edf , param_t & param )
{

  // Tan et al: (2001) https://www.ncbi.nlm.nih.gov/pubmed/11459695

  // Internight reliability and benchmark values for computer analyses
  // of non-rapid eye movement (NREM) and REM EEG in normal young
  // adult and elderly subjects.  Tan X1, Campbell IG, Feinberg I.
  // Clin Neurophysiol. 2001 Aug;112(8):1540-52.


  // Rapid Eye Movement Density is Reduced in the Normal Elderly Nato
  // Darchia, PhD1,2; Ian G. Campbell, PhD1; Irwin Feinberg

  // The demonstration that the reduced EOG power in the elderly is
  //  almost entirely due to a reduced incidence of EOG potentials is
  //  biologi- cally meaningful. If EOG amplitude had been reduced, it
  //  could have indicated an age-dependent degradation of the
  //  corneo-retinal potential difference, a reduction of eye-movement
  //  velocity, or both.23 These pos- sibilities are effectively ruled
  //  out by the finding that it is the incidence of eye-movement
  //  potentials that is reduced in the elderly; that is, elder- ly
  //  normal subjects produce eye movements during REM sleep at a
  //  lower rate than do young adults.  Reduced EMD in the elderly is
  //  presumably a degenerative (aging) rather than a developmental
  //  (maturational) age change.



}
