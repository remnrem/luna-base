
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


// RESOAP : iterative staging 

void suds_indiv_t::resoap_alter1( edf_t & edf , int epoch , suds_stage_t stage )
{
    
  // actual number of epochs
  
  const int ne_all = edf.timeline.num_epochs();
  
  // for SOAP, number of included epochs based on good signal data 
  // (i.e. as all epochs are included for SOAP targets, irrespective for 
  // ovserved stage being known or not)
  
  const int ne_included = obs_stage_valid.size();
  
  // nb. 'epoch' is 1-based , given by the user
  if ( epoch < 1 || epoch > ne_all ) 
    Helper::halt( "bad epoch value, outside range" );
  
  // some epochs may be skipped, e.g. due to signal outliers
  
  // valid epochs  : std::vector<int> epochs;   
  // all stages    : std::vector<suds_stage_t> obs_stage; 
  // valid stages  : std::vector<suds_stage_t> obs_stage_valid; 
  // same, but str : std::vector<std::string> y;   (send to lda_t() 
  
  // we need to update only y[i]; so the 'original' is kept in obs_stage(_valid)
    
  //
  // Update 'y' and check we have en
  //
  
  // epochs[] contains the codes of epochs actually present in the model/valid                                                                                                                      
  //  (in 0, ... encoding)
  
  bool updated = false;
  
  for (int i=0; i < epochs.size(); i++) 
    {
      
      // internal epoch number = epochs[i]
      // display epoch = edf.timeline.display_epoch( i )
      // nb. in SOAP context, with no restructuring of the EDF, the display epoch
      // will typically be +1 the internal epoch, i.e. not expecting discontinous
      // codes;  the user is expected to 
      
      // for y and obs_stage_valid
      int e0 = i;
      
      // for obs_stage
      int e1 = epochs[i];
      
      // for user disokay epoch (1-based)
      int e2 = edf.timeline.display_epoch( e1 );

      // update this single 'observed' stage
      if ( epoch == e2 ) 
	{
	  logger << "  changing epoch " << epoch << " from " << y[e0] << " to " << suds_t::str( stage ) << "\n";
	  y[e0] = suds_t::str( stage );
	  // obs_stage_valid[ e0 ] = stage;
	  // obs_stage[ e1 ] = stage;
	  updated = true;
	}
      
      // track what we have     
    }
  
  if ( ! updated ) 
    logger << "  no updates made: did not find epoch " << epoch << " (with valid signal data)\n";
  
  
}


void suds_indiv_t::resoap_pickN( edf_t & edf , int pick )
{
  // for evaluation of SOAP only: 

  // pick N each of 'labels', using the original observed stages
  if ( obs_stage_valid.size() != y.size() )
    Helper::halt( "cannot use RESOAP pick without original staging" );
  
  // first scrub
  for (int i=0; i < suds_t::cached.y.size(); i++)
    suds_t::cached.y[i] = suds_t::str( SUDS_UNKNOWN );
  
  const int nss = suds_t::labels.size();

  std::map<std::string,int> scounts;

  // N or more  versus exactly N
  bool exact = pick < 0 ;
  if ( exact ) pick = -pick;

  const int n = y.size();

  // Yates-Fisher shuffle to get a random ordering
  std::vector<int> a( n );
  CRandom::random_draw( a );
  
  std::set<std::string> done;
  for (int i=0; i<n; i++)
    {
      
      int p = a[i]; // random draw

      std::string ss = suds_t::str( obs_stage_valid[p] );
      if ( ss == "?" ) continue;
	       
      if ( exact )
	{
	  // only add if fewer than needed?
	  int c = scounts[ ss ];	  
	  if ( c < pick )
	    {
	      y[p] = ss;
	      ++scounts[ ss];	      
	    }
	}
      else
	{	    
	  y[p] = ss;
	  ++scounts[ ss ];
	}
      
      // done for this stage?
      if ( scounts[ y[p] ] == pick )
	done.insert( y[p] );

      // all done?
      if ( done.size() == nss ) break;	
    }

}


void suds_indiv_t::resoap( edf_t & edf , bool epoch_level_output )
{

  logger << "  re-SOAPing...\n";

  //
  // this impacts format of epoch-level output
  //

  suds_t::soap_mode = 2;
  
  //
  // Count "observed" stages
  //
  
  const int n = y.size();
  
  std::map<std::string,int> ycounts;
  for (int i=0; i<n; i++) ++ycounts[ y[i] ];
  

  //
  // requires at leasrt two stages w/ at least 3 observations, and has to 
  // be greater than the number of PSCs
  //
  
  const int required_n = 3;
  
  logger << "  epoch counts:";
  int s = 0;
  int t = 0;
  int tt = 0;
  std::map<std::string,int>::const_iterator yy = ycounts.begin();
  while ( yy != ycounts.end() )
    {
      logger << " " << yy->first << ":" << yy->second;     
      tt += yy->second;
      if ( yy->first != "?" && yy->second >= required_n ) 
	{
	  ++s;
	  t += yy->second;
	}
      ++yy;
    }
  logger << "\n";

  writer.value( "S" , s );
  writer.value( "OBS_N" , t ); // need at least 3 of each for 't'
  writer.value( "OBS_P" , t/(double)tt );

  
  bool okay = s >= 2 ;
  
  // for p predictors, require at least p+2 observations
  if ( ! ( t > nc+1 ) ) okay = false;
  
  if ( ! okay )
    {
      logger << "  not enough non-missing stages for LDA with " << nc << " predictors\n";
      writer.value( "FIT" , 0 );
      return;
    }
  
  
  //
  // Re-fit the LDA
  //
  
  Eigen::MatrixXd pp;
  
  int dummy = self_classify( NULL , &pp );

  if ( dummy == 0 )
    {
      logger << "  LDA model could not converge with the current stage proposal\n";
      writer.value( "FIT" , 0 );      
      return;
    }


  //
  // Model okay
  //

  writer.value( "FIT" , 1 );


  //
  // output stage probabilities 
  //
  
  const double epoch_sec = edf.timeline.epoch_length();
  
  std::vector<std::string> final_pred = suds_t::max( pp , lda_model.labels );
  
  summarize_kappa( final_pred , true );
  
  // actual number of epochs
  const int ne_all = edf.timeline.num_epochs();

  const int bad_epochs = summarize_stage_durations( pp , lda_model.labels , ne_all , epoch_sec );

  if ( epoch_level_output )
    summarize_epochs( pp , lda_model.labels , ne_all , edf );
  
}



int suds_indiv_t::resoap_update_pp( std::vector<suds_stage_t> * st , 
				    const double th , 
				    Eigen::MatrixXd & pp ) 
{
  
  const int rows = pp.rows();
  const int cols = pp.cols();
  if ( st->size() != rows ) 
    Helper::halt( "internal error in resoap_update_pp()" );

  std::vector<suds_stage_t> st2 = *st;
  
  // counts: (high-conf) stages, epochs 
  std::set<suds_stage_t> stgs;
  std::map<suds_stage_t,int> stgs2;
  int blanked = 0;
  for (int i=0; i<rows; i++)
    {
      stgs.insert( st2[i] );
      const double mx = pp.row(i).maxCoeff();
      if ( mx < th ) 
	{
	  ++blanked;
	  st2[i] = SUDS_UNKNOWN;
	}
      else
	stgs2[ st2[i] ]++;
    }
  
  const int kept = rows - blanked;
  const int nstg = stgs.size();
  const int nstg2 = stgs2.size();
  
  logger << " nstg, kep, blanked = " 
	 << nstg << " " << nstg2 << " " 
	 << kept << " " << blanked << "\n";
  
  
  //
  // Check we have sufficient number of high-confidence assignments;
  //

  // not sure we need this... if impacts only output (not done here)
  //suds_t::soap_mode = 2;
  
  bool okay = nstg == nstg2;
  
  //
  // requires at leasrt two stages w/ at least 3 observations, and has to 
  // be greater than the number of PSCs
  //
  
  // Rule 1; stages present each need X high-conf. values 
  // const int required_n = 3;
  
  // Rule 2: for p predictors, require at least p+2 observations
  // if ( ! ( t > nc+1 ) ) okay = false;
  
  if ( ! okay ) return 0;
  
  //
  // Re-fit the LDA
  //

  posteriors_t prediction;

  if ( 0 && suds_t::qda )
    {
      qda_t qda( suds_t::str( st2 ) , U );     
      qda_model = qda.fit( suds_t::flat_priors );
      if ( ! qda_model.valid ) return 0;
      prediction = posteriors_t( qda_t::predict( qda_model , U ) ) ; 
    }
  else
    {
      lda_t lda( suds_t::str( st2 ) , U );     
      lda_model = lda.fit( suds_t::flat_priors );      
      if ( ! lda_model.valid ) return 0;
      prediction = posteriors_t( lda_t::predict( lda_model , U ) ) ; 
    }
  
  //
  // get predictions
  //
  
  //
  // output stage probabilities 
  //
  
  // needed?
  //std::vector<std::string> pp2 = suds_t::str( suds_t::max( prediction.pp , lda_model.labels ) );
  //  if ( pp2.size() != st->size() ) Helper::halt( "internal error" );
  
  int nchanged = 0;
 
  //
  // Update
  //

  //
  // All done
  //
  
  return nchanged;
}

