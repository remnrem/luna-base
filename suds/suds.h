
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


#ifndef __SUDS_H__
#define __SUDS_H__

#include <string>
#include "stats/matrix.h"
#include <vector>
#include <set>
#include <map>
#include "stats/lda.h"

#include "eval.h"
#include "helper/helper.h"

struct edf_t;

struct param_t;


enum suds_stage_t
  {
   SUDS_WAKE = 0 ,
   SUDS_N1 = 1 ,
   SUDS_N2 = 2 ,
   SUDS_N3 = 3 ,
   SUDS_REM = 5 ,
   SUDS_ARTIFACT = 6 , 
   SUDS_UNKNOWN = 7 
  };


struct suds_indiv_t {

  suds_indiv_t() { } 

  suds_indiv_t( const std::string & id ) : id(id) { } 
  
  // wrapper to add a trainer 
  void add_trainer( edf_t & edf , param_t & param );

  // process either trainer or target;
  int proc( edf_t & edf , param_t & param , bool trainer = false );
  
  // write trainers to file
  void write( edf_t & edf , param_t & param ) const;

  // read trainer data from disk
  void reload( const std::string & filename , bool load_psd = false );

  // fit LDA, io.e. after reloading U
  void fit_lda();

  // make predictions given a different individuals signal data
  lda_posteriors_t predict( const suds_indiv_t & trainer );

  // add a prediction from one trainer
  void add( const std::string & id , const lda_posteriors_t & );
  
  // get KL weights across trainers
  Data::Vector<double> wgt_kl() const;

  // individual ID
  std::string id;

  // trainer or no?  (i.e. has manual staging?)
  bool trainer;
  
  // number of final epochs
  int nve;

  // number of spectral variables
  int nbins;
  
  // spectral data: only loaded for 'weight trainers'
  // not needed to be reloaded for standard trainers
  Data::Matrix<double> PSD;
  
  // SVD
  Data::Matrix<double> U;  // based on own data
  Data::Matrix<double> U_projected; // can be projected into this space
  Data::Vector<double> W;
  Data::Matrix<double> V;

  // Hjorth (mean/variance, per signal)
  Data::Vector<double> mean_h2, sd_h2;
  Data::Vector<double> mean_h3, sd_h3;
  // for targets only, keep epoch level Hjorths (epoch x signal)
  Data::Matrix<double> h2, h3;

  // LDA
  std::vector<std::string> y;
  lda_model_t model; // calculated on reload for trainers
  
  // staging
  std::vector<suds_stage_t> obs_stage;       // always all epochs
  std::vector<suds_stage_t> obs_stage_valid; // will match prd_stage
  std::vector<suds_stage_t> prd_stage;
  
  std::map<std::string,int> counts;

  //
  // retained epochs
  //

  std::vector<int> epochs;

  //
  // target predictions/staging
  //
  
  std::map<std::string,Data::Matrix<double> > target_posteriors;
  std::map<std::string,std::vector<suds_stage_t> > target_predictions;
  
  bool operator<( const suds_indiv_t & rhs ) const {
    return id < rhs.id;
  }
  
};


struct suds_t { 

  //  friend struct suds_indiv_t;
    
  static void attach_db( const std::string & , bool );
  
  static void score( edf_t & edf , param_t & param );

  static void set_options( param_t & param )
  {

    // number of PSC components
    nc = param.has( "nc" ) ? param.requires_int( "nc" ) : 4 ;

    // smoothing factor (multiple of SD)
    denoise_fac = param.has( "lambda" ) ? param.requires_dbl( "lambda" ) : 0.5 ; 

    // epoch-level outlier removal for trainers
    if ( param.has( "th" ) ) outlier_ths = param.dblvector( "th" );
    
    standardize_u = ! param.has( "unnorm" );

    use_best_guess = ! param.has( "no-best-guess" );

    // if target staging present, ignore 
    // e.g. if is is all 'UNKNOWN'

    ignore_target_priors = param.has( "ignore-prior" );
    
    // weights: take top N % (if 0 use all)
    wgt_percentile = param.has( "pct" ) ? param.requires_dbl( "pct" ) : 0 ;
    if ( wgt_percentile < 0 || wgt_percentile > 100 ) Helper::halt( "pct should be between 0 and 100" );

    // by default, requires 5 of each 5 epochs to include a trainer
    required_epoch_n = 5;
    if ( param.has( "req-epochs" ) ) required_epoch_n = param.requires_int( "req-epochs" );
    
    //
    // channels, w/ sample rates
    //
    // sig=A,B,C
    // lwr=10,10,10
    // upr=20,20,20
    // inc=0.25,0.25,0.25
    // sr=100,100,100
    
    verbose = param.has( "verbose" );

    if ( param.requires( "sig" ) == "*" ) Helper::halt( "requires sig to be set explicitly" );

    siglab = param.strvector( "sig" );
    
    ns = siglab.size();
    
    //
    // channel-specific options can be given
    //
    
    lwr.resize( ns , 0.5 );
    upr.resize( ns , 20 );
    fac.resize( ns , 1 );
    sr.resize( ns , 100 );
    
    if ( param.has( "lwr" ) )
      {
	lwr = param.dblvector( "lwr" );
	if ( lwr.size() != ns ) Helper::halt( "incorrect number of values for lwr" );
      }
    
    if ( param.has( "upr" ) )
      {
	upr = param.dblvector( "upr" );
	if ( upr.size() != ns ) Helper::halt( "incorrect number of values for upr" );
      }
    
    if ( param.has( "fac" ) )
      {
	fac = param.intvector( "fac" );
	if ( fac.size() != ns ) Helper::halt( "incorrect number of values for fac" );
      }
    
    if ( param.has( "sr" ) )
      {
	sr = param.intvector( "sr" );
	if ( sr.size() != ns ) Helper::halt( "incorrect number of values for sr" );
      }
    

    eannot_file = param.has( "eannot" ) ? param.value( "eannot" ) : "" ;

    eannot_prepend = param.has( "prefix" ) ? ( param.value( "prefix" ) + "_" ) : "" ;
    
    mat_dump_file = param.has( "mat" ) ? param.value( "mat" ) : "" ;

  }

  
  //
  // SUDS parameters, needed to be the same across all individuals
  //

  static bool verbose;

  static int nc;

  static int ns;
  
  static std::vector<std::string> siglab;

  static std::vector<double> lwr;

  static std::vector<double> upr;

  static std::vector<int> fac;

  static std::vector<int> sr;
  
  static double wgt_percentile;
  
  static double denoise_fac;

  static bool standardize_u;
  
  static bool use_best_guess;
  
  static bool ignore_target_priors;

  static int required_epoch_n;

  static std::vector<double> outlier_ths;

  // based on trainer mean/SD (averaged), per signal
  static std::vector<double> lwr_h2, upr_h2;
  static std::vector<double> lwr_h3, upr_h3;

  static std::string eannot_file;
  static std::string eannot_prepend;

  static std::string mat_dump_file;
  
private: 

  // trainer library
  static std::set<suds_indiv_t> bank;

  // weight-trainer library
  static std::set<suds_indiv_t> wbank;


public:

  //
  // Misc helpers 
  //
  
  static void make01( Data::Matrix<double> & r ) { 

    if ( r.dim2() != 5 ) Helper::halt( "internal error, maxpp()" );
    const int n = r.dim1();
    for (int i=0;i<n;i++)
      {
	int m = 0;
	double mx = r(i,0);
	for (int j=1;j<5;j++) 
	  if ( r(i,j) > mx ) { mx = r(i,j) ; m = j; } 
	for (int j=0;j<5;j++) r(i,j) = 0;
	r(i,m) = 1;
      } // next row/epoch

  }


  static double maxpp( const Data::Vector<double> & r ) { 
    if ( r.size() != 5 ) Helper::halt( "internal error, maxpp()" );
    double mx = r[0];    
    for (int j=1;j<5;j++) 
      if ( r[j] > mx ) mx = r[j] ;
    return mx;
  }

  static std::string max( const Data::Vector<double> & r ) { 
    if ( r.size() != 5 ) Helper::halt( "internal error, max()" );
    int m = 0;
    double mx = r[0];

    for (int j=1;j<5;j++) 
      if ( r[j] > mx ) { mx = r[j] ; m = j; } 
    if ( m == 0 ) return "N1";
    if ( m == 1 ) return "N2";
    if ( m == 2 ) return "N3";
    if ( m == 3 ) return "REM";
    return "W";
  }

  static int num( const std::string & ss ) {
    if ( ss == "N1" ) return -1;
    if ( ss == "N2" ) return -2;
    if ( ss == "N3" ) return -3;
    if ( ss == "REM" ) return 0;
    if ( ss == "W" ) return 1;
    return 2;
  }


  // downcast
  static std::string NRW( const std::string & ss ) { 
    if ( ss == "REM" ) return "R";
    if ( ss == "N1" || ss == "N2" || ss == "N3" ) return "NR";
    return "W";
  }

  static std::vector<std::string> NRW( const std::vector<std::string> & ss ) { 
    std::vector<std::string> s( ss.size() );
    for (int i=0;i<ss.size();i++) s[i] = NRW( ss[i] );
    return s;

  }

  static std::vector<suds_stage_t> type( const std::vector<std::string> & s )
  {
    std::vector<suds_stage_t> pp( s.size() );
    for (int i=0;i<s.size();i++) pp[i] = type( s[i] );
    return pp;
  }
  
  static std::vector<std::string> str( const std::vector<suds_stage_t> & s )
  {
    std::vector<std::string> pp( s.size() );
    for (int i=0;i<s.size();i++) pp[i] = str( s[i] );
    return pp;
  }

  static std::string str( const suds_stage_t & s )
  {
    if ( s == SUDS_WAKE ) return "W";
    if ( s == SUDS_N1 ) return "N1";
    if ( s == SUDS_N2 ) return "N2";
    if ( s == SUDS_N3 ) return "N3";
    if ( s == SUDS_REM ) return "REM";
    if ( s == SUDS_ARTIFACT ) return "BAD";
    if ( s == SUDS_UNKNOWN ) return "?";       
    return "?";
  }
  
  static suds_stage_t type( const std::string & s )
  {
    if ( s == "W" ) return SUDS_WAKE;
    if ( s == "N1" ) return SUDS_N1;
    if ( s == "N2" ) return SUDS_N2;
    if ( s == "N3" ) return SUDS_N3;
    if ( s == "REM" ) return SUDS_REM;
    if ( s == "BAD" ) return SUDS_ARTIFACT;
    if ( s == "?" ) return SUDS_UNKNOWN;
    return SUDS_UNKNOWN;
  }

  
  static std::map<std::string,std::map<std::string,int> > tabulate( const std::vector<std::string> & a , 
								    const std::vector<std::string> & b , 
								    const bool print = false );
    
};


#endif
