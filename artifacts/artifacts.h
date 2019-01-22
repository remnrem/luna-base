
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


#ifndef __ARTIFACTS_H__
#define __ARTIFACTS_H__

#include <string>

struct annot_t;
struct edf_t;
struct param_t;

void      rms_per_epoch( edf_t & , param_t & );

void      mse_per_epoch( edf_t & , param_t & );

void      lzw_per_epoch( edf_t & , param_t & );

annot_t * brunner_artifact_detection( edf_t & edf , 
				      const std::string & signal , 
				      const std::string & filename = "" );

annot_t * buckelmuller_artifact_detection( edf_t & edf , 
					   param_t & param , 
					   const std::string & signal , 
					   const double delta_threshold = 2.5 , 
					   const double beta_threshold = 2.0 , 
					   const double delat_lwr = 0.6 , 
					   const double delta_upr = 4.6 ,
					   const double beta_lwr = 40 , 
					   const double beta_upr = 60 ,
					   const std::string & filename = "" );

void    spike_signal( edf_t & edf , int s1 , int s2 , double wgt , const std::string & ns = "" );

#endif
