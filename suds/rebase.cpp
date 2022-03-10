
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

#include "suds.h"

#include <vector>
#include <map>
#include <set>
#include <iomanip>

#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"

#include "dirent.h"

#include "stats/eigen_ops.h"
#include "stats/lda.h"
#include "stats/statistics.h"
#include "miscmath/miscmath.h"
#include "miscmath/crandom.h"

#include "edf/edf.h"
#include "edf/slice.h"

#include "dsp/resample.h"
#include "fftw/fftwrap.h"
#include "dsp/mtm/mtm.h"
#include "dsp/tv.h"

extern logger_t logger;

extern writer_t writer;


//
// REBASE : change epoch len
//

void suds_indiv_t::rebase( edf_t & edf , param_t & param , double elen )
{
    
  // track ID (needed if caching for RESOAP)
  id = edf.id;

  // this impacts whether epochs w/ missing values are dropped or not  
  suds_t::soap_mode = 1;

  // ensure we do not call self_classify() from proc
  suds_t::self_classification = false;

  // cannot ignore existig staging in REBASE mode (in first run)
  suds_t::ignore_target_priors = false;

  // assume that we have manual staging ('true') 
  int n_unique_stages = proc( edf , param , true );
  
  // Perhaps no observed stages?
  if ( n_unique_stages < 2 )
    {
      logger << "  *** fewer than 2 non-missing stages for this individual, cannot complete REBASE\n";
      return;
    }
  
  // fit LDA: populates suds_indiv_t::model object
  fit_qlda();

  if ( ! lda_model.valid )
    {
      logger << "  *** not enough data/variability to fit LDA\n";
      return;
    }

  
  // save this old self  
  suds_indiv_t old_self = *this;
  
  // now change epoch size to target  
  edf.timeline.set_epoch( elen , elen , 0 ) ;
  
  // and re-estimate PSD assuming no known staging ('false')
  // (this will also calculate PSC, but we will ignore this... add option to skip that in proc() in future)  
  suds_t::ignore_target_priors = true;

  // also clear this , as 'summarize_epochs() will try to use it otherwise in output)
  obs_stage.clear();

  // re-process file
  n_unique_stages = proc( edf , param , true ); 

  // true means has staging (i.e. not a 'target' in the SUDS sense, but 
  // but the suds_t::ignore_target_priors means this is ignored (i.e. we do 
  // not try to reference the staging (which presumably no longer matches the epoch 
  // duration)
  
  // now project & predict into self's prior PSC space;  i.e. use same model, but will just be
  // based on PSD estimated from differently-sized epochs

  posteriors_t new_staging = predict( old_self , suds_t::qda );
  

  //
  // output stage probabilities ( new_staging.pp ) 
  //

  const double epoch_sec = edf.timeline.epoch_length();

  const int ne_all = edf.timeline.num_epochs();

  std::vector<std::string> final_pred = suds_t::max( new_staging.pp , lda_model.labels );

  const int bad_epochs = summarize_stage_durations( new_staging.pp , lda_model.labels , ne_all , epoch_sec );

  summarize_epochs( new_staging.pp , lda_model.labels , ne_all , edf );
  
}

