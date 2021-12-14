
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
// PLACE : place existing stages
//

void suds_indiv_t::place( edf_t & edf , param_t & param , const std::string & stagefile  )
{

  //
  // Get .eannot style staging data  (i.e. not expected to match EDF in duration)
  //
  
  if ( stagefile == "" ) Helper::halt( "no stages=<file> given" );
  
  if ( ! Helper::fileExists( Helper::expand( stagefile ) ) )
    Helper::halt( "problem opening " + stagefile );
  
  std::ifstream IN1( Helper::expand( stagefile ).c_str() , std::ios::in );

  std::vector<std::string> allstages;
  while ( ! IN1.eof() )
    {
      std::string line;
      Helper::safe_getline( IN1 , line );
      if ( IN1.eof() ) break;
      if ( line == "" ) continue;
      if ( line[0] == '%' ) continue;
      std::vector<std::string> tok = Helper::parse( line );
      if ( tok.size() != 1 ) Helper::halt( "expecting one stage per line" );
      allstages.push_back( line );
    }
  IN1.close();
  
  const int nstages = allstages.size();

  const int nedf = edf.timeline.num_epochs();

  logger << "  read " << nstages << " epochs from " << stagefile << "\n";
  logger << "  based on EDF, there are " << nedf << " " << edf.timeline.epoch_length() << "-s epochs\n";
  
  if ( nstages == nedf ) 
    {
      logger << "  nothing to do, epoch and EDF epoch counts are equal\n"; 
      return;
    }
  
  
  //
  // Run initial 
  //

  // // track ID (needed for caching)
  // id = edf.id;

  // // this impacts whether epochs w/ missing values are dropped or not  
  // // 0 SUDS
  // // 1 SOAP
  // // 2 RESOAP/PLACE .. i.e. allow missing
  // suds_t::soap_mode = 2;

  // // ensure we do not call self_classify() from proc
  // suds_t::self_classification = false;

  // // assume that we do /not/ have manual staging initially ('false') 
  // int n_unique_stages = proc( edf , param , false );
  

  // //
  // // Shift... start from leftmost and go all way to rightmost
  // //   stages might be longer than EDF, or shorter

  // //  STAGES           01234567
  // //     EDF           0123 
  
  // //                   0123
  // //  STAGES   0123456 7...
  // //  STAGES    012345 67..
  // //  STAGES     01234 567.
  // //  STAGES      0123 4567
  // //  STAGES       012 3456 7
  // //  STAGES        01 2345 67
  // //  STAGES         0 1234 567
  // //  STAGES           0123 4567
  // //  STAGES           .012 34567
  // //  STAGES           ..01 234567
  // //  STAGES           ...0 1234567


  // //  STAGES           01
  // //     EDF           0123 

  // //  STAGES         0 1...
  // //  STAGES           01..
  // //  STAGES           .01.
  // //  STAGES           ..01
  // //  STAGES           ...0 1

  // int estart = - ( nstages - 1 );
  // int estop   = nedf - 1 ; // inclusive

  // for (int s=estart; s<= estop; s++)
  //   {
  //     int p = estart;

  //     std::vector<std::string> trial( nedf , "." );

  //     for (int i=0; i<nstages; i++)
  // 	{
  // 	  if ( p >= 0 && p < nedf )
  // 	    {
  // 	      trial[p] = allstages[i];
  // 	      std::cout << " edf " << p << " --> stage " << i << "\n";
  // 	    }
  // 	  ++p;
  // 	}
      
  //     // now evaliate this set of stages in trial

  //   }

  // //
  // // fit LDA, and extract posteriors ( --> pp ) 
  // //

  // Eigen::MatrixXd pp;
  
  // int dummy = self_classify( NULL , &pp );
  
  // if ( dummy == 0 ) 
  //   {
  //     logger << "  *** not enough data/variability to fit LDA\n";
  //     return;
  //   }

  
  // //
  // // output stage probabilities 
  // //

  // const double epoch_sec = edf.timeline.epoch_length();

  // const int ne_all = edf.timeline.num_epochs();

  // std::vector<std::string> final_pred = suds_t::max( pp , lda_model.labels );

  // summarize_kappa( final_pred , true );

  // const int bad_epochs = summarize_stage_durations( pp , lda_model.labels , ne_all , epoch_sec );
  
  // if ( epoch_level_output )
  //   summarize_epochs( pp , lda_model.labels , ne_all , edf );



  // //
  // // RESOAP... 
  // //

  // if ( suds_t::cached.id != edf.id )
  //   Helper::halt( "need to SOAP w/ 'save' option before running RESOAP" );

  // // check that this same individual has been cached by 
  // // a previous SOAP run
  // if ( suds_t::cached.id != edf.id ) 
  //   Helper::halt( "need to SOAP w/ 'save' option before running RESOAP" );

  // // need to reset only y[]
  // // keep obs_stage[] and obs_stage_valid[] as is (i.e. any 'original' true staging)

  // //
  // // scrub all stages?
  // //
  
  // if ( param.has( "scrub" ) )
  //   {
  //     for (int i=0; i < suds_t::cached.y.size(); i++)
  // 	suds_t::cached.y[i] = suds_t::str( SUDS_UNKNOWN );            
  //     return;
  //   }
  
  // //
  // // pick N of each epoch at random?
  // //
  
  // if ( param.has( "pick" ) )
  //   {
  //     int n = param.requires_int( "pick" );
  //     suds_t::cached.resoap_pickN( edf , n );
  //     suds_t::cached.resoap( edf , param.has( "verbose" ) );
  //     return;
  //   }

  // //
  // // else, alter a single epoch
  // //
  
  // // which epoch is being updated...
  // int epoch = param.requires_int( "epoch" );
  // // ...to which stage?
  // suds_stage_t stage = suds_t::type( param.requires( "stage" ) );

  // // update and refit model based on set PSC 
  // suds_t::cached.resoap_alter1( edf , epoch , stage );
  // suds_t::cached.resoap( edf , param.has( "verbose" ) );



}

