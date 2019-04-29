
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

#ifndef __LUNA_LWPREP_H__
#define __LUNA_LWPREP_H__



struct edf_t;
struct param_t;
struct retval_t;
struct sstore_t;

#include <map>
#include <set>
#include <string>

struct lw_prep_t { 

  lw_prep_t( edf_t & edf , param_t & param );
  

  // helper functions

  
  // hard-coded routines to extract key info from the retval_t's
  
  void insert_epoch2stage( retval_t & retval , const std::string & indiv , sstore_t * ss );

  void insert_stage_summary( retval_t & retval , const std::string & indiv , sstore_t * ss );

  std::set<std::string> get_annots( retval_t & retval , const std::string & indiv );

  void insert_annot2ints( retval_t & retval , const std::string & indiv , const std::string & annot , sstore_t * ss );

  void insert_psd_band( retval_t & retval , const std::string & indiv , sstore_t * ss );

  void insert_psd_spec( retval_t & retval , const std::string & indiv , sstore_t * ss );
  
  void insert_exe_clusters( retval_t & retval , const std::string & indiv , sstore_t * ss );
  
  // for time-course power
  bool denoise;
  double lambda;
  
};

#endif

