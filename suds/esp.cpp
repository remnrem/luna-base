
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


void suds_t::read_elapsed_stages( const std::string & f )
{
  // alraady attached?
  if ( ES_probs.rows() != 0 ) return;

  // expecting format: ES PP(N1) PP(N2) PP(N3) PP(R) PP(W)
  // where ES is the prior number of elapsed sleep epochs before this one (minutes)
  // and the probabilities are based on the average in this range (i.e. up to the next ES)

  if ( ! Helper::fileExists( f ) )
    Helper::halt( "could not find ES model file " + f );
  
  std::vector<double> pp_n1, pp_n2, pp_n3, pp_r, pp_w;

  ES_mins.clear();
  
  std::ifstream IN1( f.c_str() , std::ios::in );
  
  while ( ! IN1.eof() )
    {
      std::string line;
      Helper::safe_getline( IN1 , line );
      if ( IN1.eof() ) break;
      if ( line == "" ) continue;
      if ( line[0] == '#' || line[0] == '%' ) continue;
      std::vector<std::string> tok = Helper::parse( line , "\t " );
      if ( tok.size() != 6 ) Helper::halt( "bad format for " + f );
      if ( tok[0] == "ES" ) continue;
      
      double c1,c2,c3,c4,c5,c6;
      if ( ! Helper::str2dbl( tok[0] , &c1 ) ) Helper::halt( "bad value in " + f );
      if ( ! Helper::str2dbl( tok[1] , &c2 ) ) Helper::halt( "bad value in " + f );
      if ( ! Helper::str2dbl( tok[2] , &c3 ) ) Helper::halt( "bad value in " + f );
      if ( ! Helper::str2dbl( tok[3] , &c4 ) ) Helper::halt( "bad value in " + f );
      if ( ! Helper::str2dbl( tok[4] , &c5 ) ) Helper::halt( "bad value in " + f );
      if ( ! Helper::str2dbl( tok[5] , &c6 ) ) Helper::halt( "bad value in " + f );
      if ( c1 < 0 )  Helper::halt( "bad value in " + f );
      if ( c2 < 0 || c2 > 1 ) Helper::halt( "bad value in " + f );
      if ( c3 < 0 || c3 > 1 ) Helper::halt( "bad value in " + f );
      if ( c4 < 0 || c4 > 1 ) Helper::halt( "bad value in " + f );
      if ( c5 < 0 || c5 > 1 ) Helper::halt( "bad value in " + f );
      if ( c6 < 0 || c6 > 1 ) Helper::halt( "bad value in " + f );
      
      ES_mins.push_back( c1 );
      pp_n1.push_back( c2 );
      pp_n2.push_back( c3 );
      pp_n3.push_back( c4 );
      pp_r.push_back( c5 );
      pp_w.push_back( c6 );
    }

  if ( ES_mins.size() < 1 ) Helper::halt( "could not read data from " + f );
  
  IN1.close();
  
  // P(ES|stage) should sum to 1.0 but just in case... 
  const int nbins = pp_n1.size();
  double s1 = 0 ,s2 = 0 ,s3 = 0 ,sr = 0 ,sw = 0;
  for (int i=0; i<nbins; i++)
    {
      s1 += pp_n1[i];
      s2 += pp_n2[i];
      s3 += pp_n3[i];
      sr += pp_r[i];
      sw += pp_w[i];
    }
  if ( s1 <= 0 || s2 <= 0 || s3 <= 0 || sr <= 0 || sw <= 0 ) Helper::halt( "bad format in " + f );
  
  for (int i=0; i<nbins; i++)
    {
      pp_n1[i] /= s1;
      pp_n2[i] /= s2;
      pp_n3[i] /= s3;
      pp_r[i] /= sr;
      pp_w[i] /= sw;
    }

  ES_probs = Eigen::MatrixXd::Zero( nbins , 5 );

  for (int i=0; i<nbins; i++)
    {
      ES_probs(i,0) = pp_n1[i];
      ES_probs(i,1) = pp_n2[i];
      ES_probs(i,2) = pp_n3[i];
      ES_probs(i,3) = pp_r[i];
      ES_probs(i,4) = pp_w[i];    
    }
  
  logger << "  read " << nbins << "-bin ES model from " << f << "\n";
  
}


Eigen::MatrixXd suds_t::apply_es_model( const Eigen::MatrixXd & pp ,
					const std::vector<std::string> & stg )
{
  // current posterior probailities
  Eigen::MatrixXd revised = pp ; 

  // current best-guess stage = stg (i.e. to calculate elapsed sleep) 

  // nb. if there are very large gaps in the valid record (i.e. big chunks of bad data)
  // then the elapsed sleep estimates will be off (obviously), so probably should
  // give a note that es-model=X might not be wanted in that scenario

  const int nr = pp.rows();
  
  std::vector<int> es_bin( nr );

  const int nbins = ES_mins.size();
  
  double elapsed_sleep = 0 ;

  // Note: **assumes** 30 second epochs and 5-class classification here
  
  double epoch_duration_mins = 0.5;
  
  // nb: this **assumes** that elapsed sleep bins should start at 0 
  int curr_bin = 0;
  
  for (int i=0; i<nr; i++)
    {
      
      if ( curr_bin < nbins - 1 && elapsed_sleep >= ES_mins[ curr_bin + 1 ] )
	++curr_bin;
      
      revised(i,0) *= revised(i,0) * ES_probs(curr_bin,0);
      revised(i,1) *= revised(i,1) * ES_probs(curr_bin,1);
      revised(i,2) *= revised(i,2) * ES_probs(curr_bin,2);
      revised(i,3) *= revised(i,3) * ES_probs(curr_bin,3);
      revised(i,4) *= revised(i,4) * ES_probs(curr_bin,4);

      // scale to sum to 1.0
      const double row_sum = revised(i,0) + revised(i,1) + revised(i,2) + revised(i,3) + revised(i,4);

      revised(i,0) /= row_sum;
      revised(i,1) /= row_sum;
      revised(i,2) /= row_sum;
      revised(i,3) /= row_sum;
      revised(i,4) /= row_sum;

      // get next ES value for next epoch
      if ( stg[i] != "W" ) elapsed_sleep += epoch_duration_mins;
    }

  // all done
			  
  return revised;
}
