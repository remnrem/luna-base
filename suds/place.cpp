
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

      // enforce N1, N2, N3, R, W, ?, L
      if ( line == "L" ) line = "?";
      if ( line != "?" && line != "N1" && line != "N2" && line != "N3" && line != "R" && line != "W" )
	Helper::halt( "stages=<file> lines can only be one of: N1, N2, N3, R, W, L or ?" );
      
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

  //
  // Required extent of overlap (by default, 10%)
  //
  
  double req_overlap = param.has( "overlap" ) ? param.requires_dbl( "overlap" ) : 0.1 ; 
  logger << "  requiring " << req_overlap << " proportion overlap\n\n";


  // force alignment even for equal epochs sizes?
  const bool force_align = param.has( "force" );
  
  if ( nstages == nedf && ! force_align ) 
    {
      logger << "  nothing to do, epoch and EDF epoch counts are equal\n"; 
      return;
    }

    
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

  //
  // the above may have dropped back epochs... we need to ensure that we
  // skip those when aligning the proposal epochs.. i.e if
  //   full EDF = 1000 epochs
  //   valid signals  = 950
  //   proposal = 600
  //  then we need to align the 600 against the full 1000 epochs ;;; but any 'drop' because of bad signal
  //   we will set the phenotype to ?, i.e. so that it is effectively dropped
  //
   
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
  // track valid epochs
  //

  std::map<int,int> e2e;
  for (int i=0; i< epochs.size(); i++)
    e2e[ epochs[i] ] = i ;
  
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

      //  1 2 3 4 5 6 7 8 9 10
      //  1 2   4     7 8 9
      
      
      // nedf = all EDF epochs
      // epochs.size() = all valid epochs
      // only pull valid epochs here
      
      // build full proposal:
      
      std::vector<std::string> ftrial( nedf , "?" );
      
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
	      ftrial[p] = allstages[i];       	      
       	    }
       	  ++p;
       	}
      
      const int overlap = maxe  - mine + 1;
      
      // std::cout << " FROM : " << p1 << "  ----  "
      //  		<< " ACT " << mine << " " << maxe << " OL = " << overlap << " " << overlap/(double)nedf << "\n";

      //
      // Extract only epochs that are valid 
      //

      const int ne_valid = epochs.size();
      std::vector<std::string> trial( ne_valid , "?" );
      for (int i=0; i<ne_valid; i++) trial[i] = ftrial[ epochs[i] ];
      
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
	      //std::cout << " adding " << yy->first << "\n";
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
      // track as a solution?
      //
      
      double overlap_fraction = overlap/(double)nedf;
    
      

      //
      // Model okay?
      //
      
      vec_fit.push_back( overlap_fraction >= req_overlap  ? 1 : 0  );
      
      //
      // Get alignment kappa
      //
      
      std::vector<std::string> prd = suds_t::max( pp , lda_model.labels );
      
      double kappa = -1 , kappa3 = -1;
      
      if ( overlap_fraction >= req_overlap )
	{
	  
	  kappa = MiscMath::kappa( prd ,
				   trial ,
				   suds_t::str( SUDS_UNKNOWN ) );
	  
	  kappa3 = MiscMath::kappa( suds_t::NRW( prd ) ,
				    suds_t::NRW( trial ) ,
				    suds_t::str( SUDS_UNKNOWN ) );
	  
	  matched = true;
	  
	  if ( kappa > max_kappa )
	    {
	      best_offset = p1;
	      max_kappa = kappa;
	    }	  
	  
	}
      
      // record 
      
      vec_k.push_back( kappa );
      vec_k3.push_back( kappa3 );
      
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

  int best_idx = -1;
  
  for (int i=0; i<vec_k.size(); i++)
    {
      writer.level( vec_offset[i] , "OFFSET" );

      // take first best match
      if ( matched && best_idx == -1 && vec_offset[i] == best_offset )
	best_idx = i;
      
      writer.value( "FIT" , vec_fit[i] );
      writer.value( "NS"  , vec_nst[i] );
      writer.value( "NE"  , vec_ne[i] );
      writer.value( "SS"  , vec_ss[i] );      
      
      writer.value( "OLAP_N" , vec_overlap[i] );
      writer.value( "OLAP_EDF" , vec_overlap[i]/(double)nedf );
      writer.value( "OLAP_INP" , vec_overlap[i]/(double)nstages );

      if ( vec_k[i] >= 0 )
	{
	  writer.value( "K" , vec_k[i] );
	  writer.value( "S" , vec_k[i] / max_k );
	}
      
      if ( vec_k3[i] >= 0 )
	{
	  writer.value( "K3" , vec_k3[i] );
	  writer.value( "S3" , vec_k3[i] / max_k3 );
	}
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

  logger << "\n  optimal epoch offset = " << ( best_offset >= 0 ? "+" : "-" ) << best_offset << " epochs (kappa = " << max_kappa << ")\n"
	 << "  which spans " << vec_overlap[ best_idx ] << " epochs (of " << nedf << " in the EDF, and of " << nstages << " in the input stages)\n";
  
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

