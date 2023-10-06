
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


extern logger_t logger;

extern writer_t writer;


//
// REBASE : change epoch len
//

void suds_indiv_t::rebase( edf_t & edf , param_t & param , double elen )
{

  // input, signals and original epoch length

  // analysis window
  //   length wlen
  //   shifts forward wlen - wolap
  //   ignores that windows that span > 1 stage annotation
  
  // final output - arbitrary window size
  //   based on averaging all above analysis windows in that interval


  //
  // build model using existing EPOCH settings (allowing for overlaps)
  //

  double elen0 = edf.timeline.epoch_length();
  double elap0 = edf.timeline.epoch_inc();
  
  logger << "  fitting SOAP model with epoch size " << elen0
	 << "s, and overlap of " << elap0 << "s\n";
  
  // track ID (needed if caching for RESOAP)
  id = edf.id;

  // this impacts whether epochs w/ missing values are dropped or not  
  suds_t::soap_mode = 1;

  // ensure we do not call self_classify() from proc
  suds_t::self_classification = false;

  // cannot ignore existing staging in REBASE mode (in first run)
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
  
  //
  // Get predictions 
  //
  
  posteriors_t pp = predict( *this , suds_t::qda );

  //
  // Get associated epoch times
  //

  std::vector<interval_t> etimes;
  edf.timeline.first_epoch();
  while ( 1 )
    {      
      int epoch = edf.timeline.next_epoch();      
      if ( epoch == -1 ) break;      
      interval_t interval = edf.timeline.epoch( epoch );
      etimes.push_back( interval );      
    }
  
  //
  // Get new implied epoch times
  //
  
  // now change epoch size to target  
  edf.timeline.set_epoch( elen , elen );

  std::vector<interval_t> newtimes;
  edf.timeline.first_epoch();
  while ( 1 )
    {
      int epoch = edf.timeline.next_epoch();
      if ( epoch == -1 ) break;
      interval_t interval = edf.timeline.epoch( epoch );
      newtimes.push_back( interval );
    }

  const int ne_old = etimes.size();
  const int ne_new = newtimes.size();
  
  posteriors_t qq;
  qq.pp = Eigen::MatrixXd::Zero( ne_new , suds_t::n_stages );
  // qq.cl.resize( new_new );
  // qq.cli.resize( new_new );

  //
  // Populate qq.pp 
  //

  int curr = 0; 

  for (int e=0; e<ne_new; e++)
    {
      // this new epoch
      const interval_t & interval = newtimes[e];

      // find overlapping old epochs, and calculate a weight for each
      // based on fractional overlap
      std::map<int,double> overlaps;

      // start search based on 'curr'
      // move backwards until before start of this epoch
      while ( 1 )
	{
	  if ( curr == ne_old ) --curr;
	  if ( curr == 0 ) break;
	  if ( etimes[curr].stop <= interval.start ) break;
	  --curr;
	}

      // we are now guaranteed to be before this test interval
      // move forwards
      while ( 1 )
	{
	  if ( curr == ne_old ) break;
	  if ( etimes[curr].start >= interval.stop ) break;
	  double p = interval.prop_overlap( etimes[curr] );
	  if ( p > 0 ) overlaps[ curr ] = p ;
	  ++curr;	  
	}
      
      Eigen::ArrayXd pnew = Eigen::VectorXd::Zero( suds_t::n_stages );
      
      std::map<int,double>::const_iterator oo = overlaps.begin();
      while ( oo != overlaps.end() )
	{
	  pnew += oo->second * pp.pp.row( oo->first ).array();
	  ++oo;
	}

      // scale to 1.0 
      pnew /= pnew.sum();

      // save
      qq.pp.row( e ) = pnew;
      
    }
  
  
  //
  // Outputs
  //

  const double epoch_sec = edf.timeline.epoch_length();
  
  const int ne_all = edf.timeline.num_epochs();

  std::vector<std::string> final_pred = suds_t::max( qq.pp , lda_model.labels );

  const int bad_epochs = summarize_stage_durations( qq.pp , lda_model.labels , ne_all , epoch_sec );

  summarize_epochs( qq.pp , lda_model.labels , ne_all , edf );


  // //
  // // OLD
  // //

  // // save this old self  
  // suds_indiv_t old_self = *this;
  
  // // now change epoch size to target  
  // edf.timeline.set_epoch( elen , elen );

  // // analysis will happen at window size elen
  // // but final stages will be output at (smoothed) resolution based on elen-eoverlap
  
  // // and re-estimate PSD assuming no known staging ('false')
  // // (this will also calculate PSC, but we will ignore this... add option to skip that in proc() in future)  
  // suds_t::ignore_target_priors = true;

  // // also clear this , as 'summarize_epochs() will try to use it otherwise in output)
  // obs_stage.clear();

  // // re-process file
  // n_unique_stages = proc( edf , param , true ); 

  // // true means has staging (i.e. not a 'target' in the SUDS sense, but 
  // // but the suds_t::ignore_target_priors means this is ignored (i.e. we do 
  // // not try to reference the staging (which presumably no longer matches the epoch 
  // // duration)
  
  // // now project & predict into self's prior PSC space;  i.e. use same model, but will just be
  // // based on PSD estimated from differently-sized epochs

  // posteriors_t new_staging = predict( old_self , suds_t::qda );
  

  // //
  // // output stage probabilities ( new_staging.pp ) 
  // //

  // const double epoch_sec = edf.timeline.epoch_length();

  // const int ne_all = edf.timeline.num_epochs();

  // std::vector<std::string> final_pred = suds_t::max( new_staging.pp , lda_model.labels );

  // const int bad_epochs = summarize_stage_durations( new_staging.pp , lda_model.labels , ne_all , epoch_sec );

  // summarize_epochs( new_staging.pp , lda_model.labels , ne_all , edf );
  
}

