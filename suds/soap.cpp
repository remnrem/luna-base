
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
// SOAP : Single Observation Accuracies & Probabilities
//

void suds_indiv_t::evaluate( edf_t & edf , param_t & param )
{
  
  // track ID (needed if caching for RESOAP)
  id = edf.id;

  // this impacts whether epochs w/ missing values are dropped or not  
  suds_t::soap_mode = 1;

  // ensure we do not call self_classify() from proc
  suds_t::self_classification = false;

  // verbose output
  bool epoch_level_output = param.has( "epoch" );

  // assume that we have manual staging ('true') 
  int n_unique_stages = proc( edf , param , true );
  
  //
  // Cache for RESOAP?
  //

  if ( suds_t::cache_target ) 
    {
      logger << "\n  caching " << id << " for a subsequent RESOAP\n";
      suds_t::cached = *this;
    }


  //
  // No observed stages?
  //

  if ( n_unique_stages < 2 )
    {
      logger << "  *** fewer than 2 non-missing stages for this individual, cannot complete SOAP\n";
      return;
    }


  //
  // fit LDA, and extract posteriors ( --> pp ) 
  //

  Eigen::MatrixXd pp;
  
  int dummy = self_classify( NULL , &pp );

  if ( dummy == 0 ) 
    {
      logger << "  *** not enough data/variability to fit LDA\n";
      return;
    }

  //
  // dump predictor matrix?
  //

  // to output stream
  if ( param.has( "feature-matrix" ) )
    {
      dump_predictor_matrix( edf , "" );
    }

  // as file
  if ( param.has( "dump-features" ) )
    {
      dump_predictor_matrix( edf , param.value( "dump-features" ) );
    }

  //
  // dump associations w/ stages
  //

  if ( param.has( "dump-stage-assocs" ) )
    {
      logger << "  dumping feature/SVD component stage associations to " << param.value( "dump-stage-assocs" )  << "\n";
      dump_stage_associations( param.value( "dump-stage-assocs" ) );
    }

  
  //
  // dump components
  //

  if ( param.has( "dump-svd" ) )
    {
      logger << "  dumping SVD components to " << param.value( "dump-svd" )  << "\n";
      dump_svd( param.value( "dump-svd" ) );
    }
  
  //
  // output stage probabilities 
  //

  logger << "\n";
  
  const double epoch_sec = edf.timeline.epoch_length();

  const int ne_all = edf.timeline.num_epochs();

  const std::vector<std::string> & labels = suds_t::qda ? qda_model.labels : lda_model.labels ; 
  
  std::vector<std::string> final_pred = suds_t::max( pp , labels );

  summarize_kappa( final_pred , true );

  const int bad_epochs = summarize_stage_durations( pp , labels , ne_all , epoch_sec );
  
  if ( epoch_level_output )
    summarize_epochs( pp , labels , ne_all , edf );


  //
  // Output annotations (of discordant epochs)
  //

  if ( param.has( "annot" ) )
    {
      const std::string annot_folder = param.has("annot-dir") ? param.value( "annot-dir" ) : "./";      
      write_annots( annot_folder , param.value( "annot" ) , pp , labels , ne_all , edf );
    }
  
}




int suds_indiv_t::self_classify( std::vector<bool> * included , Eigen::MatrixXd * pp )
{

  if ( ! trainer )
    Helper::halt( "can only self-classify trainers (those w/ observed staging" );
  
  // assume putative 'y' and 'U' will have been constructed, and 'nve' set
  // i.e. this will be called after proc(), or from near the end of proc() 

  //
  // fit the LDA to self
  //


  fit_qlda();

  if ( suds_t::qda && ! qda_model.valid )
    return 0;
  
  if ( (!suds_t::qda) && ! lda_model.valid )
    return 0;

  
  //
  // get predictions
  //

  posteriors_t prediction;
  if ( suds_t::qda )
    prediction = posteriors_t( qda_t::predict( qda_model , U ) ) ; 
  else
    prediction = posteriors_t( lda_t::predict( lda_model , U ) ) ; 

  // save posteriors?
  if ( pp != NULL ) *pp = prediction.pp ;


  //
  // In SOAP mode, all done (we only needed the PP)
  //
  
  if ( suds_t::soap_mode || included == NULL )
    return 1;  // SOAP only cares about a non-zero return value

  //
  // Get kappa 
  //

  double kappa = MiscMath::kappa( prediction.cl , y , suds_t::str( SUDS_UNKNOWN )  );

  included->resize( nve , false );

  
  //
  // Optionally, ask whether trainer passes self-classification kappa threshold.  If not
  // make all epochs 'bad', i.e. so that this trainer will not be used
  //
  
  if ( suds_t::self_classification_kappa <= 1 )
    {
      if ( kappa < suds_t::self_classification_kappa )
	{
	  logger << "  trainer does not meet SOAP kappa " << kappa << " < " << suds_t::self_classification_kappa << "\n";
	  return 0;  // all 'included' false at this point
	}      
    }

  
  //
  // Determine 'bad' epochs
  //

  int okay = 0;

  // hard calls

  if ( suds_t::self_classification_prob > 1 )
    {
      for (int i=0;i<nve;i++)
	{
	  (*included)[i] = prediction.cl[i] == y[i] ; 
	  if ( (*included)[i] ) ++okay;
	}
    }
  else
    {
      logger << "  using threshold of PP > " << suds_t::self_classification_prob << "\n";

      // map labels to slots in PP matrix (this might be non-standard, e.g. no REM) for a
      // given trainer, and so we cannot assume canonical slot positions
      
      std::vector<std::string> labels = suds_t::qda ? qda_model.labels : lda_model.labels;
      
      std::map<std::string,int> label2slot;
      for (int j=0;j< labels.size();j++)
	label2slot[ labels[j] ] = j ;
      
      // check PP for the observated stage
      for (int i=0;i<nve;i++)
	{
	  std::map<std::string,int>::const_iterator ii = label2slot.find( y[i] );
	  if ( ii == label2slot.end() )
	    Helper::halt( "internal error in suds_indiv_t::self_classify() , unrecognized label" );
	  
	  //	  std::cout << " i " << i << "/" << nve << "  " << prediction.pp( i , ii->second ) << " " << suds_t::self_classification_prob << "\n";
	  
	  // must match hard call::: 
	  (*included)[i] = prediction.cl[i] == y[i] ;
	  
	  if ( (*included)[i] ) 
	    ++okay;
	  else // but also, must be above threshold
	    {
	      if ( prediction.pp( i , ii->second ) >= suds_t::self_classification_prob )
		{
		  (*included)[i] = true;
		  ++okay;
		}
	      else
		(*included)[i] = false;
	    }
	}
    }

  return okay;
}

