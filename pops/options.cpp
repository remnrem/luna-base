
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

#include "pops/options.h"

#include "stats/lgbm.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "stats/eigen_ops.h"
#include "db/db.h"

extern logger_t logger;
extern writer_t writer;

#include "eval.h"

bool pops_opt_t::verbose;
double pops_opt_t::spectral_resolution;
int pops_opt_t::trim_wake_epochs;
int pops_opt_t::n_stages;
bool pops_opt_t::welch_median;

double pops_opt_t::lwr;
double pops_opt_t::upr;

std::vector<double> pops_opt_t::slope_range{ 30.0 , 45.0 } ;
double pops_opt_t::slope_th  = 3;
double pops_opt_t::slope_epoch_th = 5;

void pops_opt_t::set_options( param_t & param )
{
  
}

#endif
