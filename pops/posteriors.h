
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

#ifndef __LUNA_POPS_POSTERIORS_H__
#define __LUNA_POPS_POSTERIORS_H__

#include "stats/Eigen/Dense"

#include <string>
#include <vector>

struct edf_t;

namespace pops_posteriors {

bool add_edf_channels( edf_t * pedf ,
                       const Eigen::MatrixXd & P ,
                       const std::vector<int> & E ,
                       const std::vector<bool> * flagged ,
                       int ne_total ,
                       bool three_state ,
                       const std::string & prefix ,
                       const std::string & source_tag );

bool add_edf_signals( edf_t * pedf ,
                      const Eigen::MatrixXd & P ,
                      const std::vector<int> & E ,
                      const std::vector<bool> * flagged ,
                      int ne_total ,
                      bool three_state ,
                      double row_duration_sec ,
                      const std::string & prefix ,
                      const std::string & source_tag );

}

#endif
#endif
