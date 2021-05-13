
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

#include "stats/cpt.h"
#include "eval.h"

#include "helper/logger.h"
#include "db/db.h"

extern writer_t writer;
extern logger_t logger;

cpt_t::cpt_t( param_t & param )
{
  
  //
  // Similar input structure as for PSC; currently, 
  //   - specify a single IV 
  //   - F, CH or CH1+CH2 are the valid stratifiers currently
  //   - an ID-linked covariate table that includes a DV
  //   - a model specification (linear / logistic) 
  //   - adjacent frequencies defined in the obvious manner
  //   - adjacent spatially defined based on clocs, based on a distance threshold
  //

  //
  // Adjacency definition for EEG channel neighbours; likewise for frequency domains
  //

  double spatial_threshold = param.has( "th-spatial" ) ? param.requires_dbl( "th-spatial" ) : 0 ;
  
  double freq_threshold = param.has( "th-freq" ) ? param.requires_dbl( "th-freq" ) : 0 ; 

  //
  // Linear or logstic DV
  //

  bool dv_continuous = param.has( "linear" ) ;

  // 'value' of DV as a string ( e.g. "CASE" or "1" ) 
  std::string dv_binary = param.has( "logistic" ) ? param.value( "logistic" ) : "" ; 

  if ( ( !  dv_continuous  ) && dv_binary == "" ) 
    Helper::halt( "need to specify either 'linear' or 'logistic=DV' arguments" );
  
  
  //
  // Covariate and DV : assume a single file
  //
  
  std::string vartable = param.requires( "file" );
  
  std::string dv_name = param.requires( "dv" );
  
  //
  // Sleep metrics : assume a single file 
  //  (may be multiple files concatenated, but fine, and we can skip rows that are 'ID' ) 
  //
  
  std::string ivartable = param.requires( "iv-file" );

  std::string iv_name = param.requires( "iv" );

  
  //
  // IV outlier removal: assume that DV and covariates have already been assessed here
  //

  double outlier_threshold = param.has( "th" ) ? param.requires_dbl( "th" ) : 0 ; 

  
  //
  // Number of replicates 
  //
  
  const int nreps = param.requires_int( "nreps" );


  //
  // Channel locations
  //

  std::string clocs_file = param.has( "clocs" ) ? param.value( "clocs" ) : "" ;
  
  

  //
  // Attach sleep metrics, and define adjacent points 
  //

  

  //
  // Attach covariates, phenotype
  //
  
  
  //
  // Initial analysis of the observed data 
  //

  
  //
  // Set up permutations
  //
  

  //
  // Pemrute
  //

  
  //
  // Report outcomes
  //
  


  //
  // All done
  //


}
  

