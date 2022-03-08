
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

  //
  // Output best fit .eannot file , i.e. that matches the original EDF
  //

  const std::string ostages = param.has( "out" ) ? param.value( "out" ) : "" ; 
  
  
  //
  // get # EDF epochs (and set, if needed)
  //

  const int nedf = edf.timeline.ensure_epoched();

  logger << "  read " << nstages << " epochs from " << stagefile << "\n";
  logger << "  based on EDF, there are " << nedf << " " << edf.timeline.epoch_length() << "-s epochs\n";
  
  if ( nstages == nedf ) 
    {
      logger << "  nothing to do, epoch and EDF epoch counts are equal\n"; 
      return;
    }

  //
  // Required extent of overlap (by default, 10%)
  //
  
  double req_overlap = param.has( "overlap" ) ? param.requires_dbl( "overlap" ) : 0.1 ; 
  
  
  //
  // Initial EDF processing
  //
  
  // track ID (needed for caching)
  id = edf.id;
  
  // this impacts whether epochs w/ missing values are dropped or not  
  // 0 SUDS
  // 1 SOAP
  // 2 RESOAP/PLACE .. i.e. allow missing
  suds_t::soap_mode = 1;
  
  // // ensure we do not call self_classify() from proc
  suds_t::self_classification = false;
  
  // true = 'is a trainer' 
  int n_unique_stages = proc( edf , param , true );
  
  std::cout << "n_uniq = " << n_unique_stages << "\n";
  
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
  
  int estart = - ( nstages - 1 );
  int estop   = nedf - 1 ; // inclusive

  //
  // Track outputs
  //

  std::vector<int> vec_fit, vec_nst, vec_ne, vec_overlap, vec_offset;
  std::vector<double> vec_k, vec_k3;
  std::vector<std::string> vec_ss;
  
  double max_kappa = -1 ;
  int best_offset = -1;
  bool matched = false;
  
  //
  // Iterate over alignments
  //
  
  for (int s1=estart; s1<= estop; s1++)
    {

      // putative offset

      int p = s1;

      //
      // construct trail staging 
      //

      std::vector<std::string> trial( nedf , "?" );
      
      int cnt = 0;
      int mine = 9999999;
      int maxe = -9999999;
      int p1 = p;
      
      for (int i=0; i<nstages; i++)
       	{
       	  if ( p >= 0 && p < nedf )
       	    {
	      ++cnt;	      
	      if ( p < mine ) mine = p;
	      if ( p > maxe ) maxe = p;	      
	      trial[p] = allstages[i];       	      
       	    }
       	  ++p;
       	}
      
      const int overlap = maxe  - mine + 1;
      
      std::cout << " FROM : " << p1 << "  ----  "
		<< " ACT " << mine << " " << maxe << " OL = " << overlap << " " << overlap/(double)nedf << "\n";

      //
      // Evaluate this set of stages in trial
      //

      y = trial;
      
      const int n = y.size();     
      std::map<std::string,int> ycounts;
      for (int i=0; i<n; i++) ++ycounts[ y[i] ];
        
      //
      // requires at least two stages w/ at least 3 observations, and has to 
      // be greater than the number of PSCs
      //
      
      const int required_n = 3;
  
      int s = 0;
      int t = 0;
      int tt = 0;
      std::stringstream ss;
      std::map<std::string,int>::const_iterator yy = ycounts.begin();
      while ( yy != ycounts.end() )
	{
	  ss << ( yy != ycounts.begin() ? "," : "" ) << yy->first << ":" << yy->second;     
	  tt += yy->second;
	  if ( yy->first != "?" && yy->second >= required_n ) 
	    {
	      ++s;
	      t += yy->second;
	    }
	  ++yy;
	}

      
      //
      // Outputs 
      //

      vec_offset.push_back( p1 );
      vec_nst.push_back( s );
      vec_ne.push_back( t );	  
      vec_ss.push_back( ss.str() );
      vec_overlap.push_back( overlap );
      
      bool okay = s >= 2 ;
      
      // for p predictors, require at least p+2 observations
      if ( ! ( t > nc+1 ) )
	okay = false;
  
      if ( ! okay )
	{	  
	  vec_fit.push_back( 0 );
	  vec_k.push_back( -1 );
	  vec_k3.push_back( -1 );
			 
	  // shift to the next window
	  continue;
	}
      
      //
      // Re-fit the LDA
      //

      Eigen::MatrixXd pp;
  
      int dummy = self_classify( NULL , &pp );

      if ( dummy == 0 )
	{

	  vec_fit.push_back( 0 );	
	  vec_k.push_back( -1 );
	  vec_k3.push_back( -1 );
	  
	  continue;
	}
      
      //
      // Model okay
      //
      
      vec_fit.push_back( 1 );

      //
      // Get alignment kappa
      //
      
      std::vector<std::string> prd = suds_t::max( pp , lda_model.labels );
      
      double kappa = MiscMath::kappa( prd ,
				      trial ,
				      suds_t::str( SUDS_UNKNOWN ) );
      
      double kappa3 = MiscMath::kappa( suds_t::NRW( prd ) ,
				       suds_t::NRW( trial ) ,
				       suds_t::str( SUDS_UNKNOWN ) );

      
      vec_k.push_back( kappa );
      vec_k3.push_back( kappa3 );


      //
      // track as a solution?
      //

      double overlap_fraction = overlap/(double)nedf;

      if ( overlap_fraction >= req_overlap )
	{
	  matched = true;
	  if ( kappa > max_kappa )
	    {
	      best_offset = p1;
	      max_kappa = kappa;
	    }
	}
      
      
    }


  //
  // Outputs
  //

  double max_k = 0, max_k3 = 0;
  for (int i=0; i<vec_k.size(); i++)
    {
      if ( vec_k[i] > max_k ) max_k = vec_k[i];
      if ( vec_k3[i] > max_k3 ) max_k3 = vec_k3[i];
    }

  for (int i=0; i<vec_k.size(); i++)
    {
      writer.level( vec_offset[i] , "OFFSET" );

      writer.value( "FIT" , vec_fit[i] );
      writer.value( "NS"  , vec_nst[i] );
      writer.value( "NE"  , vec_ne[i] );
      writer.value( "SS"  , vec_ss[i] );      
      
      writer.value( "OLAP_N" , vec_overlap[i] );
      writer.value( "OLAP_P" , vec_overlap[i]/(double)nedf );

      writer.value( "K" , vec_k[i] );
      writer.value( "K3" , vec_k3[i] );

      writer.value( "S" , vec_k[i] / max_k );
      writer.value( "S3" , vec_k3[i] / max_k3 );
    }

  writer.unlevel( "OFFSET" );
  
      
  //
  // did we find an optimal point?
  //
  
  if ( ! matched )
    {
      logger << "  not able to find an optimal alignment that satisfies the overlap requirement\n" ; 
      return;
    }

  //
  // write datafiles out
  //

  if ( ostages != "" )
    {

      std::vector<std::string> trial( nedf , "?" );      
      int cnt = 0;
      int p = best_offset;
      for (int i=0; i<nstages; i++)
       	{
       	  if ( p >= 0 && p < nedf )
       	    {
	      ++cnt;	      
	      trial[p] = allstages[i];       	      
       	    }
       	  ++p;
       	}

      logger << "  writing aligned stage file (.eannot) to " << ostages << "\n";
      std::ofstream OUT1( Helper::expand( ostages ).c_str() , std::ios::out );
      for (int i=0; i<trial.size(); i++)
	OUT1 << trial[i] << "\n";      
      OUT1.close();

      
    }
  
  

  
}

