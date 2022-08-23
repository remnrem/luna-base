
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

#ifdef HAS_LGBM

#include "pops/indiv.h"

#include "helper/helper.h"
#include "helper/logger.h"
#include "stats/eigen_ops.h"
#include "db/db.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "dsp/resample.h"

#include "fftw/fftwrap.h"
#include "pdc/pdc.h"
#include <cmath>


// model:  P( elapsed sleep, prior NREM | stage )
//    elapsed sleep = duration of elapsed sleep (in 20-min bins)
//    prior NREM    = minutes of prior NREM, allowing for up to X epochs of non-NREM (10-mins bins)

Eigen::MatrixXd pops_t::ES_probs;            // P( ES, prior-NREM | stage )
std::vector<double> pops_t::ES_mins;         // total mins elapsed sleep
std::vector<double> pops_t::ES_prior_nrem; 
std::map<int,std::map<int,int> > pops_t::ES_rowmap;

// for target, either count only most likely stage, versus use weights
bool pops_t::ES_fractional_count = false;    

// as we update an epoch, use the newly-updated counts when looking at the next epochs
bool pops_t::ES_rolling = false;             

// amount of non-NREM allowed when counting prior NREM 
double pops_t::non_NREM_mins = 5;

//   N W W W N N N R N N ? ....
//                       X
//   6       5 4 3   2 1 
//     4 3 2       1
//             
//  if allowing 2 mins of non-NREM, then prior NREM = 5 epochs worth
//  if allowing 5 mins, then prior NREM = 6, etc
//  if allowing 0, prior NREM = 2


// original ES priors:
//   based on P( elapsed sleep | stage )
//    - derived from observed training data
//    - applied to targets to update posteriors, based on most likely current stage

// changes:

//   added ES_prior_nrem as well as ES_mins
//    - model P( elapsed sleep, prior NREM duration (< X non-NREM ) | stage )

//        es-weighted  --> ES_fractional_count 
//    - calculated elapsed sleep in target based on most likely stage (rather than using weights)

//        es-rolling   --> ES_rolling
//    - did not update PP derived elapsed sleep after weighting -- i.e. if updating epoch e, this should impact e+1


// i.e. prior
// ES   P(ES|N1)  P(ES|N2)   ...
// 0    0.10       0.10      ...
// 20   0.08       0.22      ...
// 40   0.02       0.35      ...
// ...
//     (probs sum to 1.0) 


//  Now we have more rows, to jointly describe a given epoch by ES but also % recent REM and NREM
//  ES   NREM   P( ES , P_NR , P_REM | N1 )   P( ES , P_NR , P_REM | N2 ) ...
//  0    0        x                             y

//  20   0     
//  20   10    
//  20   20    

//  40   0 
//  40   10 
//  40   20
//  40   30
//  40   40 
// etc


void pops_indiv_t::apply_espriors( const std::string & f )
{

  //
  // only load priors only once
  //
  
  if ( pops_t::ES_probs.rows() == 0 )
    {
      
      const std::string filename = Helper::expand( f );

      if ( ! Helper::fileExists( filename ) )
	Helper::halt( "could not open " + filename );
      
      // expecting format: ES PP(N1) PP(N2) PP(N3) PP(R) PP(W)
      // where ES is the prior number of elapsed sleep epochs before this one (minutes)
      // and the probabilities are based on the average in this range (i.e. up to the next ES)
      
      std::vector<double> pp_n1, pp_n2, pp_n3, pp_r, pp_w;
      
      pops_t::ES_mins.clear();
      
      std::ifstream IN1( filename.c_str() , std::ios::in );

      int row = 0;
      
      while ( ! IN1.eof() )
	{
	  std::string line;
	  Helper::safe_getline( IN1 , line );
	  if ( IN1.eof() ) break;
	  if ( line == "" ) continue;
	  if ( line[0] == '#' || line[0] == '%' ) continue;
	  std::vector<std::string> tok = Helper::parse( line , "\t " );
	  if ( tok.size() != 7 ) Helper::halt( "bad format for " + filename );
	  if ( tok[0] == "ES" ) continue;
	  
	  double c1,c2,c3,c4,c5,c6,c7;
	  if ( ! Helper::str2dbl( tok[0] , &c1 ) ) Helper::halt( "bad value in " + filename );
	  if ( ! Helper::str2dbl( tok[1] , &c2 ) ) Helper::halt( "bad value in " + filename );
	  if ( ! Helper::str2dbl( tok[2] , &c3 ) ) Helper::halt( "bad value in " + filename );
	  if ( ! Helper::str2dbl( tok[3] , &c4 ) ) Helper::halt( "bad value in " + filename );
	  if ( ! Helper::str2dbl( tok[4] , &c5 ) ) Helper::halt( "bad value in " + filename );
	  if ( ! Helper::str2dbl( tok[5] , &c6 ) ) Helper::halt( "bad value in " + filename );
	  if ( ! Helper::str2dbl( tok[6] , &c7 ) ) Helper::halt( "bad value in " + filename );
	  
	  if ( c1 < 0 )  Helper::halt( "bad value in " + filename );
	  if ( c2 < 0 )  Helper::halt( "bad value in " + filename );
	  if ( c3 < 0 || c3 > 1 ) Helper::halt( "bad value in " + filename );
	  if ( c4 < 0 || c4 > 1 ) Helper::halt( "bad value in " + filename );
	  if ( c5 < 0 || c5 > 1 ) Helper::halt( "bad value in " + filename );
	  if ( c6 < 0 || c6 > 1 ) Helper::halt( "bad value in " + filename );
	  if ( c7 < 0 || c7 > 1 ) Helper::halt( "bad value in " + filename );


	  // track definition
	  pops_t::ES_mins.push_back( c1 );
	  pops_t::ES_prior_nrem.push_back( c2 );
	  std::cout << " c1, c2 = " << c1 <<" " << c2 << " --> " << row << "\n";
	  pops_t::ES_rowmap[ (int)c1 ][ (int)c2 ] = row;
	  ++row;
	  
	  // read in probs
	  pp_n1.push_back( c3 );
	  pp_n2.push_back( c4 );
	  pp_n3.push_back( c5 );
	  pp_r.push_back( c6 );
	  pp_w.push_back( c7 );
	}
      
      if ( pops_t::ES_mins.size() < 1 )
	Helper::halt( "could not read data from " + filename );
      
      IN1.close();


      //
      // we can assume that P(ES,prior NREM|stage) should sum to 1.0 but just in case...
      //
      
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

      if ( s1 <= 0 || s2 <= 0 || s3 <= 0 || sr <= 0 || sw <= 0 )
	Helper::halt( "bad format in " + f );
      
      for (int i=0; i<nbins; i++)
	{
	  pp_n1[i] /= s1;
	  pp_n2[i] /= s2;
	  pp_n3[i] /= s3;
	  pp_r[i] /= sr;
	  pp_w[i] /= sw;
	}

      //
      // Construct prior prob map
      //
      
      pops_t::ES_probs = Eigen::MatrixXd::Zero( nbins , 5 );
      
      // Note -- align columns in same order as P for POPS:
      //   W R N1 N2 N3
      for (int i=0; i<nbins; i++)
	{
	  pops_t::ES_probs(i,0) = pp_w[i];	  
	  pops_t::ES_probs(i,1) = pp_r[i];
	  pops_t::ES_probs(i,2) = pp_n1[i];
	  pops_t::ES_probs(i,3) = pp_n2[i];
	  pops_t::ES_probs(i,4) = pp_n3[i];
	}
      
      logger << "  read " << nbins << "-bin ES model from " << filename << "\n";
      
    }
  
  
  //
  // apply ES priors model
  //

  // inputs: P  =  posteriors P( stage , prior NREM | signals ) 
  //         S  =  assigned stage

  logger << "  applying ES prior model...\n";
  

  //
  // revise the current posterior probailities, given elapsed sleep priors
  //
  
  Eigen::MatrixXd revised = P ;
  
  
  // use the current best-guess stage (S) to calculate elapsed sleep
  
  // nb. if there are very large gaps in the valid record (i.e. big chunks of bad data)  
  // then the elapsed sleep estimates will be off (obviously), so probably should
  // give a note that es-model=X might not be wanted in that scenario
  
  const int nr = revised.rows();
  
  const int nbins = pops_t::ES_mins.size();

  // initialize these -- both 0 at the start
  double elapsed_sleep = 0 ;
  double recent_nrem = 0;
  
  // Note: **assumes** 30 second epochs and 5-class classification here
  
  double epoch_duration_mins = 0.5;
  
  // nb: this **assumes** that elapsed sleep bins should start at 0 
  int curr_bin = 0;
  
  for (int i=0; i<nr; i++)
    {
      
      // get bin numbers :: nb. using hard-coding here
      
      int es_bin = floor( elapsed_sleep > 360 ? 360 : elapsed_sleep / 20.0 );
      int nrem_bin = floor( recent_nrem > 60 ? 60 : recent_nrem / 10.0 );
      
      int es_min = es_bin * 20;
      int nrem_min = nrem_bin * 10;

      if ( pops_t::ES_rowmap.find( es_min ) == pops_t::ES_rowmap.end() )
	Helper::halt( "internal error in finding ES bin(1)" );
      
      std::cout << " E = " << i << "\t" << S[i] << " -->  nrem_bin = " << es_min << " " << nrem_min << " || " << recent_nrem << "\n";

      if ( pops_t::ES_rowmap[ es_min ].find( nrem_min ) == pops_t::ES_rowmap[ es_min ].end() )
	Helper::halt( "internal error in finding NR bin(2)" );
      
      curr_bin = pops_t::ES_rowmap[ es_min ][ nrem_min ];

      // update probs
      revised(i,0) *= revised(i,0) * pops_t::ES_probs(curr_bin,0);
      revised(i,1) *= revised(i,1) * pops_t::ES_probs(curr_bin,1);
      revised(i,2) *= revised(i,2) * pops_t::ES_probs(curr_bin,2);
      revised(i,3) *= revised(i,3) * pops_t::ES_probs(curr_bin,3);
      revised(i,4) *= revised(i,4) * pops_t::ES_probs(curr_bin,4);

      // scale to sum to 1.0
      const double row_sum = revised(i,0) + revised(i,1) + revised(i,2) + revised(i,3) + revised(i,4);
      
      revised(i,0) /= row_sum;
      revised(i,1) /= row_sum;
      revised(i,2) /= row_sum;
      revised(i,3) /= row_sum;
      revised(i,4) /= row_sum;

      // get next ES value for next epoch
      if ( S[i] != pops_stage_t::POPS_WAKE )
	elapsed_sleep += epoch_duration_mins;
      
      // get next NR value for next epoch
      int j = i ; 
      bool first_NREM = false;      
      int recent_nrem_epochs = 0;      
      const double epoch_mins = 0.5 ; // nb. hard-coded
      int allowance = pops_t::non_NREM_mins / epoch_mins ;
      int nonNREM = 0;
      while ( 1 )
	{
	  if ( j == 0 ) break;
	  
	  // move back in time
	  --j;
	  
	  // is NREM?
	  const bool is_nrem = S[j] == POPS_N1 || S[j] == POPS_N2 || S[j] == POPS_N3 ;

	  // track whether we've yet hit any NREM
	  if ( is_nrem ) first_NREM = true;
	  
	  // track stage
	  if ( is_nrem )
	    ++recent_nrem_epochs;
	  else if ( first_NREM ) // only start counting once we've hit some NREM 
	    ++nonNREM;
	  
	  if ( nonNREM > allowance ) break;
	    
	}
      
      recent_nrem = recent_nrem_epochs * epoch_mins;

      
    }

  // all done  
  P = revised;
  
}


//
// stand-alone function to make es-priors file from training data
//

void pops_t::make_espriors( param_t & param )
{

  //   1) read binary data --> (optioally write as text) --> calc/write es-priors
  //   2) read text 

  //   extract :: read from 'data' (binary) --> write to 'text' 
  // ! extract :: read from prior text
  
  const bool extract = param.has( "extract" );
  
  const std::string data_file = param.has( "data" ) ? Helper::expand( param.value( "data" ) ) : "";
  
  const std::string text_file = param.has( "text" ) ? Helper::expand( param.value( "text" ) ) : "";

  if ( data_file == "" && text_file == "" )
    Helper::halt( "no input specified ('data' or 'text' args)");

  //
  // populate S, E, etc from binary file
  //

  if ( data_file != "" )
    {
      // scan data file, but only extract SS (not feature matrix)
      load1_stages_only( data_file );
      
      if ( text_file != "" )
	{
	  logger << "  writing epoch/stage info to " << text_file << "\n";

	  std::ofstream O1( text_file.c_str(), std::ios::out );
	  
	  for (int e=0; e<S.size(); e++)
	    O1 << E[e] << "\t" << S[e] << "\n";
	  
	  O1.close();
	}
      
    }
  else
    {
      //
      // read from text file
      //

      if ( ! Helper::fileExists( text_file ) )
	Helper::halt( "could not open " + text_file + " for reading" );
      
      S.clear();
      E.clear();
      
      std::ifstream IN1( text_file.c_str() , std::ios::in );
      while ( IN1 )
	{
	  int epoch;
	  int stg;
	  IN1 >> epoch >> stg;
	  if ( IN1.bad() || IN1.eof() ) break;	  
	  E.push_back( epoch );
	  S.push_back( stg );
	}
      IN1.close();
      logger <<"  read " << S.size() << " epochs\n";
    }
  
  
  //
  // calculate and write ES-priors?
  //
  
  if ( ! param.has( "es-priors" ) )
    {
      logger << "  no es-priors=<file> specified, so quitting\n";
      return;
    }
  
  const std::string espriors_file = param.value( "es-priors" );
  
  //
  // parameters 
  //

  // bin size (mins)
  const double tbin = param.has( "es-min" ) ? param.requires_dbl( "es-min" ) : 20 ;
  const double nr_tbin = param.has( "nr-min" ) ? param.requires_dbl( "nr-min" ) : 10 ;
  
  // max time (mins)
  const double tmax = param.has( "es-max" ) ? param.requires_dbl( "es-max" ) : 380 ;
  const double nr_tmax = param.has( "nr-max" ) ? param.requires_dbl( "nr-max" ) : 60 ;
  
  // intercept (i.e. to avoid 0-weight probs for any cell)
  const double c    = param.has( "es-c" )   ? param.requires_dbl( "es-c" ) : 0.01 ;

  //
  // calculate and report : just needs S and E
  //
  
  write_elapsed_sleep_priors( espriors_file , tbin, tmax , nr_tbin, nr_tmax, c );
  
  // all done
  
}


void pops_t::write_elapsed_sleep_priors( const std::string & f ,
					 double es_tbin, double es_tmax,
					 double nrem_tbin , double nrem_tmax, 
					 double c )
{
  
  // Given OBSERVED S and E, calculate overall elapsed sleep prior distribution 
  // and then P(ES|stage) from training data;  save to a file that 
  // POPS es-priors=X can use during testing 

  // c    : constant, i.e. to ensure some weight all values
  // tbin : bin size (e.g. 20 mins default) 
  // tmax : maximum limit (default = 400 mins, 6.6 hrs)

  //
  
  const int ne = S.size();

  // new indiv can be inferred if E[i] <= E[i-1] 
  
  int prior_epoch = 999999;

  const double epoch_mins = pops_opt_t::epoch_inc / 60.0 ;

  double elpased_sleep_mins;
  
  std::map<int,std::map<int,std::map<int, double> > > ES; // stg -> ES-bin -> NREM-bin -> count 
  
  // last bin is that value plus (e.g. 400+)
  const int es_nbins = floor( es_tmax / es_tbin ) + 1;

  const int nrem_nbins = floor( nrem_tmax / nrem_tbin ) + 1;
  
  for (int i=0; i<ne; i++)
    {
      
      // new indiv? reset all counters
      if ( E[i] < prior_epoch )
	elpased_sleep_mins = 0;
      
      // track elapsed sleep
      elpased_sleep_mins += S[i] == POPS_WAKE ? 0 : epoch_mins ;
      
      // 'recent NREM duration'
      int recent_nrem_epochs = 0;

      // number of epochs of non-NREM to allow
      int allowance = pops_t::non_NREM_mins / epoch_mins ; 
      int nonNREM = 0;

      int j = i ; 

      bool first_NREM = false;
      
      while ( 1 )
	{
	  if ( j == 0 ) break;
	  
	  // move back in time
	  --j;
	  
	  // is NREM?
	  const bool is_nrem = S[j] == POPS_N1 || S[j] == POPS_N2 || S[j] == POPS_N3 ;

	  // track whether we've yet hit any NREM
	  if ( is_nrem ) first_NREM = true;
	  
	  // track stage
	  if ( is_nrem )
	    ++recent_nrem_epochs;
	  else if ( first_NREM ) // only start counting once we've hit some NREM 
	    ++nonNREM;
	  
	  // have we gone over?
	  if ( nonNREM > allowance )
	    break;			  
	}
      
      double recent_nrem_mins = recent_nrem_epochs * epoch_mins;
      
      // W, R, N1, N2, N3
      if ( S[i] >= 0 && S[i] <= 5 ) 
	{      
	  
	  //
	  // elapsed sleep
	  //
	  
	  int es_b = floor( elpased_sleep_mins / es_tbin ) ;

	  int nrem_b = floor( recent_nrem_mins / nrem_tbin );

	  if ( es_b < es_nbins && nrem_b < nrem_nbins )
	    {
	      // record
	      ES[ S[i] ][ es_b ][ nrem_b ]++;
	    }     
	}
      
      // track epoch for next 
      prior_epoch = E[i];
      
    }


  //
  // total number of bins  
  //

  // not all bins are possible: count the possible ones
  int tot_bins =0; 
  for (int es_bin = 0 ; es_bin < es_nbins ; es_bin++ )
    for (int nrem_bin = 0 ; nrem_bin < nrem_nbins ; nrem_bin++ )
      {
        const double es_minutes = es_bin * es_tbin;
        const double nrem_minutes = nrem_bin * nrem_tbin;
        if ( nrem_minutes <= es_minutes ) ++tot_bins;
      }
  
  
  Eigen::MatrixXd P = Eigen::MatrixXd::Zero( tot_bins , pops_opt_t::n_stages );

  int row = 0;
  
  for (int es_bin = 0 ; es_bin < es_nbins ; es_bin++ )
    for (int nrem_bin = 0 ; nrem_bin < nrem_nbins ; nrem_bin++ )
      {

	// skip 'impossible' options
	const double es_minutes = es_bin * es_tbin;
	const double nrem_minutes = nrem_bin * nrem_tbin;
	if ( nrem_minutes > es_minutes ) continue;
	
	// else, add as a row
	for (int ss = 0; ss<5; ss++)
	  P( row , ss ) = ES[ ss ][ es_bin ][ nrem_bin ];
	
	// next tot-bin
	++row;
	//std::cout << " count row = " << row << " of "<< tot_bins << "\n";
      }
  
  
  // normalize within stage
  for (int ss = 0; ss<5; ss++)
    P.col( ss ) /= P.col( ss ).sum();
  
  // add offset, re-normalize
  for (int ss = 0; ss<5; ss++)
    {      
      P.col( ss ).array() += c;
      P.col( ss ) /= P.col( ss ).sum();
    }
  
  
  //
  // Write to a file
  //
  
  std::string filename = Helper::expand( f );

  std::ofstream OUT1( filename.c_str() , std::ios::out );
  
  logger << "  writing P( elapsed sleep | stg ) to " << filename << "\n";
  
  OUT1 << "ES\t"
       << "RECENT_NR\t"
       << "PP(N1)\t"
       << "PP(N2)\t"
       << "PP(N3)\t"
       << "PP(R)\t"
       << "PP(W)\n";
  
  // to get output format: N1 N2 N3 R W
  std::vector<int> sidx = { 2 , 3 , 4 , 1 , 0 };
  std::map<int,std::vector<double> > nn;
  
  row = 0;
  for (int es_bin = 0 ; es_bin < es_nbins ; es_bin++ )
    for (int nrem_bin = 0 ; nrem_bin < nrem_nbins ; nrem_bin++ )
      {	
	// skip 'impossible' rows
	if ( nrem_bin * nrem_tbin > es_bin * es_tbin ) continue;

	OUT1 << es_bin * es_tbin << "\t"
	     << nrem_bin * nrem_tbin ;
	
	for (int s=0; s<5; s++)
	  OUT1 << "\t" << P(row,sidx[s]);

	nn[ nrem_bin ].resize( 5 , 0 );
	for (int s=0; s<5; s++)
	  nn[ nrem_bin ][s] += P(row,sidx[s]);

	OUT1 << "\n";
	
	++row;
      }
  
  OUT1.close();


  for (int nrem_bin = 0 ; nrem_bin < nrem_nbins ; nrem_bin++ )
    {
      std::cout << nrem_bin * nrem_tbin ;
      for (int s=0; s<5; s++)
	std::cout <<"\t" << nn[ nrem_bin ][s];
      
      std::cout << "\n";
    }
}


#endif
