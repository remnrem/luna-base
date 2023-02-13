
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

// pops_t         ::  main arguments, etc, and holds LGBM 
// pops_spec_t    ::  specifies which features, channels, etc to use (EDF --> feature matrix)
// pops_indiv_t   ::  data for one person

//  1) process trainers, one at a time, and write out 'core/raw' features (level 1) to a binary format

//  2a) can concatenate binary feature matrix files w/ a simple 'cat'
//  2b) load combined trainer (& validation dataset?) --> but how is it used?
//  2c) (within person) add any a) temporal smoothing, b) normalization, c) SVD , d) combining --> level 2 feature set
//  2d) fit and save model

//  3a) prediction - load model, create level 1 feature set
//  3b) create level 2 feature set (optionally loading SVD matrices, as needed) 
//  3c) predict


// new points
//  --> features:: as SUDS but perhaps also cross-spectral terms?
//  --> option:: level 2 feature sets can include e.g. 2 epochs:  i.e. to predict label N2-N2, or N2-N3 etc
//       ? idea - that this might help at transitions? 
//
//  ... - N2 - N3 - N3 - N2 - ...
//        [      ]
//            [       ]
//                 [       ]
//

#ifdef HAS_LGBM

#include "pops/pops.h"
#include "pops/indiv.h"
#include "pops/options.h"
#include "pops/spec.h"

#include "lgbm/lgbm.h"
#include "miscmath/miscmath.h"
#include "miscmath/crandom.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "stats/eigen_ops.h"
#include "stats/statistics.h"
#include "db/db.h"
#include "dsp/tv.h"

#include "helper/zfile.h"

extern logger_t logger;
extern writer_t writer;

//pops_opt_t pops_t::opt;
lgbm_t pops_t::lgbm;
std::string pops_t::lgbm_model_loaded = "";
pops_specs_t pops_t::specs;

std::map<std::string,double> pops_t::range_mean;
std::map<std::string,double> pops_t::range_sd;

std::vector<std::string> pops_t::labels5 = { "W" , "R" , "N1" , "N2" , "N3" };
std::vector<std::string> pops_t::labels3 = { "W" , "R" , "NR" };


//
// create a level 2 feature library and save
// i.e. for trainer and/or validation library
//  however, with the dump=training|test 
//  option set, instead of fitting the LGBM model, 
//  here we just compute the full feature matrix
//  and write to a file; either the training, validation
//  or test dataset;  if the test dataset, we will still read
//  the SVD solutions;  if training, it computes them, but dumps
//  all (i.e. will be both training + validation, but ordered
//   - as per usual, the SVD solution is based on training + validation
//     


void pops_t::make_level2_library( param_t & param )
{

  //
  // Inputs
  //

  // always require an explicit, single data file
  const std::string data_file   = param.requires( "data" );
  
  // features: add path as needed

  std::string feature_file  = ".";
  if ( param.has( "features" ) )
    feature_file = param.value( "features" );
  else if ( pops_opt_t::pops_root != "" )
    feature_file = pops_opt_t::pops_root + ".ftr";
  if ( feature_file != "." )
    feature_file = pops_t::update_filepath( feature_file );
  if ( feature_file == "." ) 
    Helper::halt( "POPS requires a feature file, via lib or features args" );
  
  std::string model_file  = ".";
  if ( param.has( "model" ) )
    model_file = param.value( "model" );
  else if ( pops_opt_t::pops_root != "" )
    model_file = pops_opt_t::pops_root + ".mod";
  if ( model_file != "." )
    model_file = pops_t::update_filepath( model_file );
  
  std::string conf_file  = ".";
  if ( param.has( "conf" ) )
    conf_file = param.value( "conf" );
  else if ( param.has( "config" ) )
    conf_file = param.value( "config" );
  else if ( pops_opt_t::pops_root != "" )
    conf_file = pops_opt_t::pops_root + ".conf";
  if ( conf_file != "." )
    conf_file = pops_t::update_filepath( conf_file );

  std::string ranges_file  = ".";
  if ( param.has( "ranges" ) )
    ranges_file = param.value( "ranges" );
  else if ( pops_opt_t::pops_root != "" && pops_opt_t::if_root_apply_ranges )
    ranges_file = pops_opt_t::pops_root + ".ranges";
  if ( ranges_file != "." )
    ranges_file = pops_t::update_filepath( ranges_file );
  
  std::string espriors_file  = ".";
  if ( param.has( "priors" ) )
    espriors_file = param.value( "priors" );
  else if ( pops_opt_t::pops_root != "" && pops_opt_t::if_root_apply_espriors )
    espriors_file = pops_opt_t::pops_root + ".priors";
  if ( espriors_file != "." )
    espriors_file = pops_t::update_filepath( espriors_file );
  

  //
  // Misc set up
  //

  bool is_trainer = true;
  
  //
  // Either dump feature matrix, or run LGBM
  //

  // mode 1: load, make lvl-2 matrix, (dump ranges), train model, save
  // mode 2: load, make lvl-2 matrix, dump to a file
  // mode 3: load, make lvl-2 matrix, dump ranges
  
  const bool run_mode2 = param.has( "dump" );
  
  std::string dump_file = "";
  
  const bool run_mode3 = param.has( "ranges-only" );

  // main use:
  const bool run_mode1 = ! ( run_mode2 || run_mode3 );

  // in normal use, we require a model file
  if ( run_mode1 )
    {
      if ( model_file == "." )
	Helper::halt( "POPS requires a model file to be specified, via lib or model" );      
    }

  // if dumping, is this for trainers or test cases?
  if ( run_mode2 )
    {            
      if ( param.value( "dump" ) == "test" )
	is_trainer = false;      
      else if ( param.value( "dump" ) != "training" )
       	Helper::halt( "'dump' must be set to either 'training' or 'test'" );
      
      dump_file = param.requires( "file" );
    }
  

  //
  // feature specifications (used to generate l1-data)
  //

  pops_t::specs.read( feature_file );
  
  
  //  
  // validation dataset for LGBM (expected in main set, but will be partitioned off)
  //
   
  if ( param.has( "hold-outs" ) )
    load_validation_ids( param.value( "hold-outs" ) );
  else if ( param.has( "validation" ) )
    load_validation_ids( param.value( "validation" ) );
  
	   
  //  
  // get previous data: assume single file, concatenated
  // this will populate: X1, S, E and Istart/Iend
  // validation IDs will be at the end
  //
  
  load1( data_file );
  

  //
  // Check for NaNs
  //
  
  pops_nan_report_t missing_vars( X1 );
  if ( missing_vars.any() )
    {
                                    
      std::vector<std::string> pp = pops_t::specs.col_label;
      std::map<int,int>::const_iterator cc = missing_vars.cols.begin();
      while ( cc != missing_vars.cols.end() )
	{
	  logger  << "  ** warning: " << pp[cc->first ] << " has " << cc->second << " missing values\n";
	  ++cc;
	}
      
      const int nt = Istart.size();
      for (int i=0; i<nt; i++)
	{
	  for (int e = Istart[i] ; e <= Iend[i] ; e++) 
	    {
	      if ( missing_vars.rows.find( e ) != missing_vars.rows.end() )
		{
		  logger << "  ** warning: " << I[ i ]  
			 << " has an epoch with " << missing_vars.rows[e] << " missing value(s)\n";		  
		}
	    }
	}
    }


  //
  // expand X1 to include space for level-2 features
  //

  X1.conservativeResize( Eigen::NoChange , pops_t::specs.na );

  
  // 
  // summarize training and validation dataset counts  
  //

  report_counts();


  //  
  // derive level-2 stats (jointly in *all* indiv. training + valid ) 
  //  -- typically trainer, but allow for test individuals here if we
  //     are using the --pops command not to train, but with dump=test
  //     (this will instruct level2() to read the SVD rather than create)
  //
  
  level2( is_trainer );


  //
  // prune trainer rows (atfer level2(), i.e. keeping SVD and smoothing in place
  // but drop these before adding other weights, etc and fitting the model
  // only prune trainer rows, not validation rows
  //

  if ( pops_opt_t::sample_fixed_n )
    {
      sample_fixed_n();
      logger << "  after pruning training epochs:\n";
      report_counts();
    }
  
  //
  // states
  //
  
  if ( pops_opt_t::run_stage_associations )
    stage_association();

  //
  // output the the feature matrix 
  //

  if ( run_mode2 )
    {
      dump_matrix( dump_file );
      return;
    }

  
  //
  // Write ranges (means/SDs) to a .range file
  //

  if ( ranges_file != "." )
    {

      // this does a) whole sample, b) per indiv
      dump_ranges( ranges_file );
      
      // all done?
      if ( run_mode3 ) // ranges-only
	return;
    }
  
  
  //
  // LGBM config (user-specified, or default for POPS)
  //
  
  if ( conf_file == "." )
    lgbm.load_pops_default_config();
  else    
    lgbm.load_config( conf_file );

  //
  // Set number of iterations?
  //
  
  if ( param.has( "iterations" ) ) 
    lgbm.n_iterations = param.requires_int( "iterations" );
  else if ( param.has( "iter" ) ) 
    lgbm.n_iterations = param.requires_int( "iter" );
  
  
  //  
  // Stage (label) weights?    
  //
  
  std::vector<double> wgts(  pops_opt_t::n_stages , 1.0 );
  
  if ( param.has( "weights" ) ) // order must be W, R, NR...
    {
      wgts = param.dblvector( "weights" );
      if ( wgts.size() != pops_opt_t::n_stages )
	Helper::halt( "expecting " + Helper::int2str( pops_opt_t::n_stages ) + " stage weights" );
      logger << "  read " << wgts.size() << " weights:";
      if      ( pops_opt_t::n_stages == 5 ) logger << " W, R, N1, N2, N3 =";
      else if ( pops_opt_t::n_stages == 3 ) logger << " W, R, NR =";
      for (int i=0; i<wgts.size(); i++) logger << " " << wgts[i] ;
      logger << "\n";
    }
  
  lgbm_label_t weights( pops_opt_t::n_stages == 5 ? pops_t::labels5 : pops_t::labels3 , wgts );

  
  //
  // Train model
  //
 
  fit_model( model_file , weights );
  

  //
  // Write elapsed-sleep priors?
  //
  
  if ( espriors_file != "." )
    write_elapsed_sleep_priors( espriors_file );      

  
}


//
// derive level 2 stats (from pops_t::specs)
//

void pops_t::level2( const bool training , const bool quiet )
{
  
  // go through all level2 blocks in order
  
  for (int s=0; s<pops_t::specs.specs.size(); s++)
    {
      pops_spec_t spec = pops_t::specs.specs[s];

      std::string l2ftr = pops_t::specs.ftr2lab[ spec.ftr ];
	    
      if ( pops_t::specs.lvl2.find( l2ftr ) == pops_t::specs.lvl2.end() )
	continue;
      
      // from cols

      std::string from_block = Helper::toupper( spec.arg[ "block" ] );
      std::string to_block = spec.block;
      const bool inplace = from_block == to_block;      
      std::vector<int> from_cols = pops_t::specs.block_cols( from_block , pops_t::specs.na );
      std::vector<int> to_cols = pops_t::specs.block_cols( to_block , pops_t::specs.na );

      if ( ! quiet ) 
	{
	  logger << "   - adding level-2 feature " << l2ftr << ": "; 
	  
	  if ( from_cols.size() > 0 )  // time track does not have "from" cols
	    logger << from_block << " (n=" << from_cols.size() << ") ";
	  
	  logger << "--> " 
		 << to_block << " (n=" << to_cols.size() << ", cols:" 
		 << to_cols[0] << "-" << to_cols[to_cols.size()-1] << ") \n";
	}
      
      // for (int j=0;j<from_cols.size();j++) logger << " " << from_cols[j] ;
      // logger << "\n";
      
      // for (int j=0;j<to_cols.size();j++) logger << " " << to_cols[j] ;
      // logger << "\n\n";

      const int nfrom = from_cols.size();
      const int nto   = to_cols.size();
      const int ne    = X1.rows();     
      const int ni    = Istart.size();
      
      //
      // SMOOTH
      //
      
      if ( spec.ftr == POPS_SMOOTH )
	{
	  const int hwin = spec.narg( "half-window" );
	  const int fwin = 1 + 2 * hwin;
	  const double mwin = 0.05;  // taper down to this weight (from 1 at center)

	  // these should always match
	  if ( nfrom != nto )
	    Helper::halt( "internal error (2) in level2()" );
	  
	  // need to go person-by-person
	  for (int j=0; j<nfrom; j++)
	    {
	      Eigen::VectorXd D = X1.col( from_cols[j] ) ;
	      for (int i=0; i<ni; i++)
		{
		  int fromi = Istart[i];
		  int sz    = Iend[i] - Istart[i] + 1;
		  D.segment( fromi , sz ) = eigen_ops::tri_moving_average( D.segment( fromi , sz ) , fwin , mwin );
		}
	      X1.col( to_cols[j] ) = D;
	    }
	}


      //
      // DENOSIE 
      //
      
      if ( spec.ftr == POPS_DENOISE )
	{
	  const double lambda = spec.narg( "lambda" );
	  
	  // these should always match
	  if ( nfrom != nto )
	    Helper::halt( "internal error (2) in level2()" );
	  
	  // need to go person-by-person
	  for (int j=0; j<nfrom; j++)
	    {
	      Eigen::VectorXd D = X1.col( from_cols[j] ) ;
	      for (int i=0; i<ni; i++)
		{
		  int fromi = Istart[i];
		  int sz    = Iend[i] - Istart[i] + 1;
		  const double sd = eigen_ops::sdev( D.segment( fromi , sz ) );
		  dsptools::TV1D_denoise( D.segment( fromi , sz ) , lambda * sd );		  
		}
	      // place back 
	      X1.col( to_cols[j] ) = D;
	    }	  
        }


      //
      // DERIV
      //
      
      if ( spec.ftr == POPS_DERIV )
        {
	  const int hw = spec.narg( "half-window" ) ;
	  
	  // these should always match
	  if ( nfrom != nto )
	    Helper::halt( "internal error (2) in level2()" );
	  
	  for (int j=0; j<nfrom; j++)
	    {
	      
	      Eigen::VectorXd D = X1.col( from_cols[j] ) ;

	      // need to go person-by-person
	      for (int i=0; i<ni; i++)
		{
		  int fromi = Istart[i];
		  int sz    = Iend[i] - Istart[i] + 1;
		  
		  // place back as slope of [ -hw , X , +hw ] epoch interval
		  eigen_ops::deriv( D.segment( fromi , sz ) , hw );
		  
		}
	      // place back 
	      X1.col( to_cols[j] ) = D;
	    }	  

	}

      //
      // CUMUL
      //
      
      if ( spec.ftr == POPS_CUMUL )
	{
	  
	  // -1 and +1  for neg and pos accumulators
	  // 0 = norm (default) .. scale to 0..1 first then add; then scale 0..1 again
	  // 2 = abs() values

	  int ctype = 0;
	  if ( spec.arg[ "type" ] == "pos" ) ctype = 1;
	  else if ( spec.arg[ "type" ] == "neg" ) ctype = -1;
	  else if ( spec.arg[ "type" ] == "abs" ) ctype = 2; // abs
	  
	  // these should always match
	  if ( nfrom != nto )
	    Helper::halt( "internal error (2) in level2()" );
	  
	  // need to go person-by-person
	  for (int j=0; j<nfrom; j++)
	    {

	      Eigen::VectorXd D = X1.col( from_cols[j] ) ;

	      for (int i=0; i<ni; i++)
		{
		  int fromi = Istart[i];
		  int sz    = Iend[i] - Istart[i] + 1;
		  
		  // place back as a 0..1 variable 
		  eigen_ops::accumulate( D.segment( fromi , sz ) , ctype );
		  
		}
	      // place back 
	      X1.col( to_cols[j] ) = D;
	    }	  
        }


      //
      // NORM (within person)
      //
            
      if ( spec.ftr == POPS_NORM )
	{

	  const double win = spec.narg( "winsor" );
	  if ( win < 0 || win > 0.5 )
	    Helper::halt( "winsor should be between 0 and 0.5" );
	  
          // these should always match                                                                                                     
          if ( nfrom != nto )
            Helper::halt( "internal error (2) in level2()" );
	  
          // need to go person-by-person
          for (int j=0; j<nfrom; j++)
            {
              Eigen::VectorXd D = X1.col( from_cols[j] ) ;
              for (int i=0; i<ni; i++)
                {
                  int fromi = Istart[i];
                  int sz    = Iend[i] - Istart[i] + 1;
		  eigen_ops::robust_scale( D.segment( fromi , sz ) , true , true , win , true );
                }
              X1.col( to_cols[j] ) = D;
            }

        }
      

      //
      // RESCALE (within person)
      //
            
      if ( spec.ftr == POPS_RESCALE )
	{

	  // by default, use 95% percentile
	  const double pct = spec.has( "pct" ) ? spec.narg( "pct" ) : 0.95;   

	  if ( pct <= 0 || pct >= 1 ) 
	    Helper::halt( "pct should be between 0 and 1" );
	  
	  // segments for estimating pct (no overlap) 	  
	  const int nsegs = spec.has( "n" ) ? spec.narg( "n" ) : 50 ;
	  
          // these should always match
	  if ( nfrom != nto )
            Helper::halt( "internal error (2) in level2()" );
	  
          // need to go person-by-person
          for (int j=0; j<nfrom; j++)
            {
              Eigen::VectorXd D = X1.col( from_cols[j] ) ;
              for (int i=0; i<ni; i++)
                {
                  int fromi = Istart[i];
                  int sz    = Iend[i] - Istart[i] + 1;
		  D.segment( fromi , sz ) = eigen_ops::percentile_scale( D.segment( fromi , sz ) , pct , nsegs );
                }
              X1.col( to_cols[j] ) = D;
            }
	  
        }

      //
      // SVD  - done differently for trainers versus targets
      //
      
      if ( spec.ftr == POPS_SVD )
	{

	  const int nc = spec.narg( "nc" );

	  // allow for 'path' option to modify where this file is                                                        
	  const std::string wvfile = pops_t::update_filepath( spec.arg[ "file" ] );
	  
	  // copy to a temporary
	  Eigen::MatrixXd D = Eigen::MatrixXd::Zero( ne , nfrom );	  
	  for (int j=0; j<nfrom; j++) D.col(j) = X1.col( from_cols[j] ) ;
	  
	  // mean-center (within each individual)
	  for (int i=0; i<ni; i++)
	    eigen_ops::scale( D.middleRows( Istart[i] , Iend[i] - Istart[i] + 1 ) , true , false );
	  
	  // trainer SVD
	  if ( training )
	    {
	      Eigen::BDCSVD<Eigen::MatrixXd> svd( D , Eigen::ComputeThinU | Eigen::ComputeThinV );
	      Eigen::MatrixXd U = svd.matrixU();
	      Eigen::MatrixXd V1 = svd.matrixV();
	      Eigen::VectorXd W1 = svd.singularValues();
	      
	      // copy U back to X1
	      for (int j=0; j<nto; j++) X1.col( to_cols[j] ) = U.col(j);
	      
	      // save W and V to a file (i.e. for use later in prediction models, to project
	      // test cases into this space)

	      if ( ! quiet ) 
		logger << "   - writing SVD W and V to " << wvfile << "\n";
	      
	      std::ofstream OUT1( wvfile.c_str() , std::ios::out );
	      OUT1 << V1.rows() << " " << nc << "\n";
	      for (int i=0;i<V1.rows(); i++)
		for (int j=0;j<nc; j++)
		  OUT1 << " " << V1(i,j) ;
	      OUT1 << "\n";
	      for (int j=0;j<nc; j++)
		OUT1 << " " << W1[j] ;
	      OUT1 << "\n";
	      OUT1.close();
	      
	    }
	  else
	    {
	      // projection of target 
	      
	      if ( V.find( wvfile ) == V.end() ) // do once
		{

		  if ( ! quiet )
		    logger << "   - reading SVD W and V from " << wvfile << "\n";
		  
		  if ( ! Helper::fileExists( wvfile ) )
		    Helper::halt( "cannot find " + wvfile + "\n (hint: add a 'path' arg to point to the .svd file" );
		  std::ifstream IN1( wvfile.c_str() , std::ios::in );
		  int nrow, ncol;
		  IN1 >> nrow >> ncol;
		  
		  Eigen::MatrixXd V0 = Eigen::MatrixXd::Zero( nrow , ncol );
		  Eigen::MatrixXd W0 = Eigen::MatrixXd::Zero( ncol , ncol );
		  if ( ncol != nc ) Helper::halt( "internal mismatch in SVD nc" );

		  for (int i=0;i<nrow; i++)
		    for (int j=0;j<ncol; j++)
		      IN1 >> V0(i,j) ;
		  for (int j=0;j<ncol; j++)
		    {
		      IN1 >> W0(j,j) ;
		      W0(j,j) = 1.0 / W0(j,j); // nb. take 1/W here
		    }
		  IN1.close();
		  
		  V[ wvfile ] = V0;
		  W[ wvfile ] = W0;
		}	    

	      //
	      // line-up
	      //

	      if ( D.cols() != V[wvfile].rows() )
		Helper::halt( "projection file does not align with number of features - please check this has not been swapped/modified\n" + wvfile );
	      
	      //
	      // project 
	      //
	      
	      Eigen::MatrixXd U_proj = D * V[ wvfile ] * W[ wvfile ];
	      
	      // copy back
	      for (int j=0; j<nto; j++) 
		X1.col( to_cols[j] ) = U_proj.col(j);

	    }
	}
      
      
      //
      // TIME 
      //
      
      if ( spec.ftr == POPS_TIME )
	{
	  
          const int order = spec.narg( "order" );

          // need to go person-by-person                                                                                                                                                    
	  for (int i=0; i<ni; i++)
	    {
	      int fromi = Istart[i];
	      int sz    = Iend[i] - Istart[i] + 1;
	      
	      X1.block( fromi , to_cols[0] , sz , nto ) = pops_t::add_time_track( sz , order );
	    }
	  
	}
      

    }
  
  

  
  //
  // Select FINAL feature set  (pops_t::specs.nf)
  //
  
  if ( pops_t::specs.nf != pops_t::specs.final2orig.size() )
    Helper::halt( "internal error (1) in level2()" );
  
  std::map<int,int>::const_iterator ff = pops_t::specs.final2orig.begin();
  while ( ff != pops_t::specs.final2orig.end() )
    {
      if ( ff->first > ff->second )
	Helper::halt( "internal error in level2()" );
      X1.col( ff->first ) = X1.col( ff->second );
      ++ff;
    }
  X1.conservativeResize( Eigen::NoChange , pops_t::specs.nf );
  

  //
  // Dump final feature matrix
  //

  if ( 0 )
    {
      std::ofstream OUT2( "X1.mat" , std::ios::out );
      OUT2 << X1 << "\n";
      OUT2.close();
    }
  
}


//
// fit and save a LGBM model
//

void pops_t::fit_model( const std::string & modelfile , 
			const lgbm_label_t & weights )
{

  // training   X1.topRows( nrows_training )   -->   S1
  // validation X1.bottomRows( nrows_validation )   --> S2
  
  // split stages 
  std::vector<int> S1 = S;
  S1.resize( nrows_training );

  std::vector<int> S2;
  for (int i=nrows_training; i<nrows_training+nrows_validation; i++)
    S2.push_back(S[i]);
  
  // training data
  lgbm.attach_training_matrix( X1.topRows( nrows_training ) );

  // labels
  lgbm.attach_training_labels( S1 );
  
  // training weights
  lgbm.add_label_weights( lgbm.training , &lgbm.training_weights, weights );
  
  // update w/ indiv-level weights?  calls lgbm.add_block_weights() 
  for (int i=0; i < pops_opt_t::iweights.size(); i++)
    attach_indiv_weights( pops_opt_t::iweights[i] , true ); // T --> training

  // apply final weights
  lgbm.apply_weights( lgbm.training , &lgbm.training_weights );
         
  // validation data?
  if ( nrows_validation ) 
    {
      lgbm.attach_validation_matrix( X1.bottomRows( nrows_validation ) );

      // labels
      lgbm.attach_validation_labels( S2 );
      
      // valdation  weights
      lgbm.add_label_weights( lgbm.validation , &lgbm.validation_weights, weights );
      
      // update w/ indiv-level weights?  calls lgbm.add_block_weights() 
      for (int i=0; i<pops_opt_t::iweights.size(); i++)
	attach_indiv_weights( pops_opt_t::iweights[i] , false ); // F --> validation
            
      // apply final weights
      lgbm.apply_weights( lgbm.validation , &lgbm.validation_weights );

    }

  // save weights?
  if ( pops_opt_t::dump_model_weights )
    dump_weights();

  // fit model
  lgbm.create_booster();
  
  // all done, write model
  lgbm.save_model( modelfile );
  
}



//
// attach a final LGBM model
//

// void pops_t::load_model( param_t & param )
// {
//   if ( ! lgbm_model_loaded != param.requires( "model" ) )
//     {
//       lgbm.load_config( param.requires( "config" ) );
//       lgbm.load_model( param.requires( "model" ) );
//       lgbm_model_loaded = param.requires( "model" );
//     }
// }

void pops_t::outliers( const Eigen::VectorXd & x ,
		       const double th ,
		       const std::vector<int> & staging , 
		       std::vector<int> * staging2 )
{


  double sum = 0;
  double sumsq = 0;
  int n = 0;
  const int n1 = x.size();
  
  // based mean/SD on staging[] good calls
  for (int i=0; i<n1; i++)
    {
      if ( staging[i] != POPS_UNKNOWN )
	{
	  sum += x[i];
	  sumsq += x[i] * x[i];
	  n++;
	}
    }
  
  if ( n < 3 ) return;
  const double mean = sum/(double)n;
  const double meansq = mean * mean ; 
  const double sd = sqrt( sumsq / (double)(n-1) - ( n/(double)(n-1) ) * meansq );
  const double upr = mean + th * sd ;
  const double lwr = mean - th * sd ;

  // update staging2
  for (int i=0; i<n1; i++)
    {    
      if ( (*staging2)[i] != POPS_UNKNOWN )
        {
	  if ( x[i] < lwr || x[i] > upr )
	    (*staging2)[i] = POPS_UNKNOWN;	  
	}
    }
  
}


void pops_t::from_single_target( const pops_indiv_t & indiv )
{
  X1 = indiv.X1;
  S = indiv.S;
  E = indiv.E;  
  Istart.resize( 1 , 0 );
  Iend.resize( 1 , S.size() - 1 );
}

void pops_t::copy_back( pops_indiv_t * indiv )
{
  indiv->X1 = X1;
}


pops_stats_t::pops_stats_t( const std::vector<int> & obs_ , 
			    const std::vector<int> & pred_ ,
			    const int nstages ,
			    const int type , 
			    const int ostage )
{

  // save either 3-class or 5-class stats
  // i.e. inputs obs and pred may be 5 or 3-class
  n = nstages;
  
  // any restrictions of epochs to look at? 
  // type:
  //   0 OAO  all epochs 
  //   1 AAA only epochs with similar flanking observed stages (i.e. 'consistent' sleep)  A-A-A
  //   2 AAX only left-epochs at a transition (i.e. if the following obs epoch is not the same)  A-B
  //   3 XAA only right-epochs at a transition (i.e. if the prior obs epoch was not the same)  B-A
  //   4 XAX only 'singleton' epochs flanked by different epochs on both sides
  //   5 TRN any transitions

  //  further, if ostage != -1, then only look at epochs with that obs stage type
  
  std::vector<int> obs; 
  std::vector<int> pred;
  
  const int ne = obs_.size();

  if ( type == 0 && ostage == -1 ) 
    {
      obs = obs_;
      pred = pred_;
    }
  else 
    {
      for (int i=0; i<ne; i++)
	{
	  const bool left_disc = i != 0 && obs_[i-1] !=obs_[i] ;
	  const bool right_disc = i < ne-1 && obs_[i+1] != obs_[i] ;
	  const bool left_right_disc = i == 0 || i == ne - 1 ? false : obs_[i-1] != obs_[i+1] ;
	  
	  bool add = true;
	  
	  if ( type == 1 ) // A-A-A
	    add = ! ( left_disc || right_disc ) ;
	  else if ( type == 2 ) // A-A-X 
	    add = right_disc && ! left_disc ;
	  else if ( type == 3 ) // X-A-A
	    add = left_disc && ! right_disc ;
	  else if ( type == 4 ) // X-A-X
	    add = left_disc && right_disc ;
	  else if ( type == 5 ) // TRN (not AAA)
	    add = left_disc || right_disc;
	  
	  // restrict to a particular class of observed stages too?
	  if ( ostage != -1 && ostage != obs_[i] ) 
	    add = false;

	  if ( add ) 
	    {
	      obs.push_back( obs_[i] );
	      pred.push_back( pred_[i] );
	    }	  
	}
    }

  // track n 
  nobs = obs.size();

  // only calculate stats if at least 10 obs of this type
  if ( nobs < 10 ) return;
      
  // other metrics: full set
  if ( ostage == -1 && type == 0 )  
    {
      // kappa
      kappa = MiscMath::kappa( obs , pred , POPS_UNKNOWN );
      
      std::vector<int> l5 = { 0 , 1 , 2 , 3 , 4 };
      std::vector<int> l3 = { 0 , 1 , 2 };
      
      acc = MiscMath::accuracy( obs , pred , 
			      POPS_UNKNOWN , 
			      n == 5 ? &l5 : &l3 ,
			      &precision , &recall , &f1 , 
			      &macro_precision , 
			      &macro_recall , 
			      &macro_f1 , 
			      &avg_weighted_precision , 
			      &avg_weighted_recall ,
			      &avg_weighted_f1 ,
			      &mcc  );
    }
  else // we only need accuracy for the restricted sets for now
    acc = MiscMath::accuracy( obs , pred , POPS_UNKNOWN );
  
  
}


std::map<int,std::map<int,int> > pops_t::tabulate( const std::vector<int> & a , 
						   const std::vector<int> & b , 
						   const bool print  )
{

  // assume : a = obs = rows
  //        : b = prd = cols
  
  std::map<int,std::map<int,int> > res;
  
  const int n = a.size();
  
  if ( n != b.size() ) 
    Helper::halt( "internal error: unequal vectors in tabulate()" );

  // includes unknown stages POPS_UNKNOWN, '?' in table
  //  (but these will be removed from kappa and other stats)

  std::set<int> uniq;
  for (int i=0;i<n;i++)
    {      
      res[ a[i] ][ b[i] ]++;
      uniq.insert( a[i] );
      uniq.insert( b[i] );
    }

  std::map<int,double> rows, cols;
  double tot = 0;
  std::set<int>::const_iterator uu = uniq.begin();
  while ( uu != uniq.end() )
    {
      std::set<int>::const_iterator jj = uniq.begin();
      while ( jj != uniq.end() )
	{
	  if ( res.find( *uu ) == res.end() )
	    res[ *uu ][ *jj ] = 0;
	  else
	    {
	      std::map<int,int> & rjj = res.find(*uu)->second;
	      if ( rjj.find( *jj ) == rjj.end() )
		res[ *uu ][ *jj ] = 0;
	    }	      
	  
	  // col/row marginals
	  rows[ *uu ] += res[ *uu ][ *jj ] ;
	  cols[ *jj ] += res[ *uu ][ *jj ] ;
	  tot += res[ *uu ][ *jj ];
	  
	  ++jj;
	}
      ++uu;
    }
  
  
  if ( print )
    {

      logger << "\t Pred:";
      std::set<int>::const_iterator uu = uniq.begin();
      while ( uu != uniq.end() )
	{
	  logger << "\t" << pops_t::label( (pops_stage_t)(*uu) );
	  ++uu;
	}
      logger << "\tTot\n";	
      
      logger << "  Obs:";
      uu = uniq.begin();
      while ( uu != uniq.end() )
	{
	  logger << "\t" <<  pops_t::label( (pops_stage_t)(*uu) ); 

	  std::set<int>::const_iterator jj = uniq.begin();
	  while ( jj != uniq.end() )
	    {
	      logger << "\t" << res[ *uu ][ *jj ];
	      ++jj;
	    }	  
	  
	  // row sums
	  logger << "\t" << Helper::pp( rows[ *uu ]/tot );
	  logger << "\n";
	  ++uu;
	}


      // col sums
      logger << "\tTot:";
      std::set<int>::const_iterator jj = uniq.begin();
      while ( jj != uniq.end() )
	{
	  logger << "\t" << Helper::pp( cols[ *jj ]/tot );
	  ++jj;
	}
      logger << "\t1.00\n\n";

      //
      // conditional probabilties  / res[][] / row[]
      //   --->  Pr ( predicted | observed )  P_COND_OBS
      //         Pr ( observed | predicted )  P_COND_PRED
      
      uu = uniq.begin();
      while ( uu != uniq.end() )
        {
	  writer.level( pops_t::label( (pops_stage_t)(*uu) ) , "OBS" );
	  std::set<int>::const_iterator jj = uniq.begin();
          while ( jj != uniq.end() )
	    {
	      writer.level( pops_t::label( (pops_stage_t)(*jj) ) , "PRED" );
	      writer.value( "N" , res[ *uu ][ *jj ] );
	      if ( rows[ *uu ] > 0 ) 
		writer.value( "P_COND_OBS" , res[ *uu ][ *jj ] / rows[ *uu ] );
	      if ( cols[ *jj ] > 0 ) 
		writer.value( "P_COND_PRED" , res[ *uu ][ *jj ] / cols[ *jj ] );
	      ++jj;
	    }
	  writer.unlevel( "PRED" );
	  ++uu;
	}
      writer.unlevel( "OBS" );
    }
  
  return res;
}


Eigen::MatrixXd pops_t::add_time_track( const int nr , const int tt )
{
  if ( nr <= 0 || tt <= 0 ) Helper::halt( "internal error in add_time_track()" );
  
  Eigen::MatrixXd T = Eigen::MatrixXd::Zero( nr , tt );
  
  for (int r=0; r<nr; r++)
    for ( int c=0; c<tt; c++)
      T(r,c) = pow( ( r / (double)nr ) - 0.5 , c+1 );
   
  return T;
  
}


void pops_t::load_validation_ids( const std::string & f )
{
  holdouts.clear();
  
  if ( ! Helper::fileExists( Helper::expand( f ) ) )
    Helper::halt( "could not open " + f );
  
  std::ifstream IN1( Helper::expand( f ).c_str() , std::ios::in );
  while ( ! IN1.eof() ) 
    {
      std::string id;
      IN1 >> id;
      if ( id == "" || IN1.eof() ) break;
      holdouts.insert( id );
    }
  IN1.close();

  logger << "  read " << holdouts.size() 
	 << " validation dataset individuals from " 
	 << f << "\n";
  
}

std::string pops_t::update_filepath( const std::string & f )
{
  
  if ( f == "" ) Helper::halt( "empty file name" );
  std::string f2 = Helper::expand( f );

  if ( pops_opt_t::pops_path == "" ) return f2;
  
  // add a path before hand (unless we're already given an absolute path)
  if ( f2[0] != globals::folder_delimiter )
    f2 = Helper::expand( pops_opt_t::pops_path + globals::folder_delimiter + f2 );

  return f2;
}


void pops_t::dump_matrix( const std::string & f )
{
  std::string dfile = Helper::expand( f );
  logger << "  dumping feature matrix to " << dfile << "\n";
  
  gzofstream Z1( dfile.c_str() , std::ios_base::out );
  
  Z1 << "SS";
  std::vector<std::string> labels = pops_t::specs.select_labels();
  for (int i=0; i< labels.size(); i++) 
    Z1 << "\t" << labels[i];       
  Z1 << "\n";
  
  for (int i=0; i<X1.rows(); i++)
    {
      // epoch label (i.e. stage)
      Z1 << pops_t::label( (pops_stage_t)S[i] );

      // epoch features
      for (int j=0; j< X1.cols(); j++)
	Z1 << "\t" << X1(i,j);
      Z1 << "\n";
    }
  Z1.close();
}


void pops_t::read_ranges( const std::string & f )
{

  if ( ! Helper::fileExists( f ) )
    Helper::halt( "could not open " + f );

  // expecting an exact line up w/ feature file and ranges
  // but we do not test this explicitly - - i.e. as we might read
  // this prior to building the feature set
  
  std::ifstream IN1( f.c_str() , std::ios::in );
  // header
  std::string str0, str1, str2, str3;
  IN1 >> str0 >> str1 >> str2 >> str3;
  if ( str0 != "ID" || str1 != "VAR" || str2 != "MEAN" || str3 != "SD" )
    Helper::halt( "bad format for " + f + "\n -- expecting columns ID, MEAN and SD" );
  
  while ( 1 )
    {
      std::string id;
      std::string varname;
      std::string str_mean , str_sd;

      double mean, sd;
      IN1 >> id >> varname >> str_mean >> str_sd;
      if ( IN1.eof() || IN1.bad() ) break;
      // only read initial overall values
      if ( id != "." ) break;
      if ( varname == "" ) continue;
      
      // skip if invalid ranges/means
      if ( str_mean == "." || Helper::iequals( str_mean, "nan" ) || Helper::iequals( str_mean, "-nan" ) )
	continue;
      
      if ( str_sd == "." || Helper::iequals( str_sd, "nan" ) || Helper::iequals( str_sd, "-nan" ) || str_sd == "0" )
	continue;

      if ( ! Helper::str2dbl( str_mean , &mean ) ) continue;
      if ( ! Helper::str2dbl( str_sd , &sd ) ) continue;
      
      if ( sd < 1e-6 ) continue; 
      
      range_mean[ varname ] = mean ;
      range_sd[ varname ] = sd ;
    }
  logger << "  read " << range_mean.size() << " valid feature mean/SD ranges from " << f << "\n";
  IN1.close();
}


void pops_t::dump_ranges( const std::string & f )
{

  // dump ranges a) overall (ID == ".")
  // and then person by person (trainer)
  // this will be called after a level2 library
  // construction

  std::string rfile = Helper::expand( f );
  
  logger << "  dumping ranges to " << rfile << "\n";
  
  std::ofstream O1( rfile.c_str() , std::ios_base::out );
  
  O1 << "ID\tVAR\tMEAN\tSD\n";

  const int nrow = X1.rows();
  const int ncol = X1.cols();

  std::vector<std::string> labels = pops_t::specs.select_labels();

  for (int i=0; i<ncol; i++)
    {
      double mean = X1.col(i).mean();
      double sd = sqrt((X1.col(i).array() - mean ).square().sum()/(nrow-1));
      O1 << ".\t"
	 << labels[i] << "\t"
	 << mean << "\t"
	 << sd << "\n";

      // if ( labels[i] == "HJORTH.EMG.V1" ) 
      // 	{
      // 	  std::cout <<"  MEAN = " << mean << " " << sd << "\n";
      // 	  std::cout <<" Ne = " << X1.rows() << "\n";
      // 	  std::cout << X1.col(i) << "\n\n";
      // 	}

    }

  //
  // now indiv-by-indiv
  //
  
  const int nt = Istart.size();

  //  std::cout << " stts " << nrow <<" " << ncol << " " << nt << "\n";
  
  for (int i=0; i<nt; i++)
    {
      
      // pull out data for this trainer only
      int fromi = Istart[i];
      int sz    = Iend[i] - Istart[i] + 1;
      Eigen::MatrixXd XI = X1.block( fromi , 0 , sz, ncol ); 
      
      // repeat as above (by col)
      for (int j=0; j<ncol; j++)
	{	  
	  //std::cout << "i , j = " << i << "\t" << j << "\n";
	  double mean = XI.col(j).mean();
	  double sd = sqrt((XI.col(j).array() - mean ).square().sum()/(sz-1));
	  
	  O1 << I[i] << "\t"
	     << labels[j] << "\t"
	     << mean << "\t"
	     << sd << "\n";
	}
      
      // next trainer
    }

  O1.close();

}

void pops_t::report_counts()
{
  std::map<int,int> ss_training, ss_valid;
  for (int i=0; i<nrows_training; i++) 
    ss_training[S[i]]++;
  for (int i=nrows_training; i<nrows_validation+nrows_training; i++) 
    ss_valid[S[i]]++;
  
  logger << "  nT=" << Istart.size() - ni_validation << " training individuals, "
	 << "nV=" << ni_validation << " (of " << holdouts.size() << " listed) validation individuals\n";
  logger << "  stage epoch counts:\n";
  std::map<int,int>::const_iterator ii = ss_training.begin();
  while ( ii != ss_training.end() )
    {
      logger << "  " << pops_t::label( (pops_stage_t)ii->first ) << "\t" 
	     << " train = " << ii->second << "\t"
	     << " validation = " << ss_valid[ ii->first ] << "\n";
      ++ii;
    }

}

void pops_t::sample_fixed_n()
{  

  // select epochs from trainers to keep
  std::map<int,std::vector<int> > all;
  for (int i=0; i<nrows_training; i++)
    all[ S[i] ].push_back( i );
  
  // ensure we have enough of each as requested
  for (int i=0; i<pops_opt_t::fixed_n.size(); i++)
    if ( all[i].size() < pops_opt_t::fixed_n[i] )
      {
	logger << "  *** requested " << pops_opt_t::fixed_n[i] << " "
	       << pops_t::label( (pops_stage_t)i ) << " epochs, but only observed " << all[i].size() << "\n" ;
	Helper::halt( "stopping - change fix=W,R,N1,N2,N3 options and re-run" );
      }

  // pick the requested subset 
  std::set<int> inc;
  int new_trainer_rows = 0;
  for (int i=0; i<pops_opt_t::fixed_n.size(); i++)
    {
      const int nobs = all[i].size();
      const int nreq = pops_opt_t::fixed_n[i];
      std::set<int> selected;
      while ( 1 )
	{
	  if ( selected.size() == nreq ) break;
	  const int idx = CRandom::rand( nobs ) ;
	  const int r = all[i][idx];
	  selected.insert( r );
	  inc.insert( r );
	}
      //      logger << "  got " << selected.size() << " stage " << i << "\n";
      new_trainer_rows += selected.size();
    }
  
  // make new versions of S, X1, E
  // also I/Iend
  // keep nrows_validation as neeeded
  // total n epoch = nrows_training + nrows_validation (stacked top/bottom)

  // std::vector<int> S  stages   
  // X1                  features   
  // E                   epoch number
  // nrows_training	   
       
       // I Iend              startstop of each indivi
       // keep I/Iend abd I as is - i.e. all trainer indivs
       // retained, even if some contribute 0 epochs Iend[i] = 
				   // weights not yet constructed so nothing to do there  
							 //       // advance IID if at last epoch of current indiv?


  // copy over

  std::vector<int> Sx = S;
  std::vector<int> Ex = E;
  Eigen::MatrixXd X1x = X1;
  std::vector<std::string> Ix = I;
  std::vector<int> Istartx = Istart;
  std::vector<int> Iendx = Iend;
  
  // wipe all 
  X1 = Eigen::MatrixXd::Zero( new_trainer_rows + nrows_validation , X1x.cols() );
  S.clear();
  S.resize( new_trainer_rows + nrows_validation );
  E.clear();
  E.resize( new_trainer_rows + nrows_validation );

  // copy features/stages/epochs
  int idx = 0;
  std::set<int>::const_iterator ii = inc.begin();
  while ( ii != inc.end() )
    {
      X1.row(idx) = X1x.row( *ii );
      S[idx] = Sx[ *ii ];
      E[idx] = Ex[ *ii ];
      ++idx;
      ++ii;
    }

  // copy indiv-IDs
  I.clear();
  Istart.clear();
  Iend.clear();

  // old number of indivs
  const int nix = Ix.size();

  int iid_idx = 0;
  int iid_start = -1;
  int iid_end = 0;

  int new_idx = 0;
  
  for (int i=0; i<nrows_training; i++)
    {
      
      // this row included?      
	if ( inc.find( i ) != inc.end() )
	  {
	    // only first row for this new person
	    if ( iid_start == -1 ) iid_start = new_idx;
	    iid_end = new_idx;
	    ++new_idx;
	  }      
      
      // advance to IID if at last epoch of current indiv?
      if ( i == Iendx[ iid_idx ] )
	{

	  // anything to add?
	  if ( iid_start != -1 )
	    {
	      Istart.push_back( iid_start );
	      Iend.push_back( iid_end );
	      I.push_back( Ix[iid_idx] );
	      // reset
	      iid_start = -1;
	    }

	  // advance to next person 
	  if ( iid_idx < nix - 1 ) 
	    ++iid_idx;
	}
    }

  // add validation data back in
  X1.block( new_trainer_rows , 0 , nrows_validation , X1.cols() ) = X1x.block( nrows_training , 0 , nrows_validation , X1.cols() );

  for (int i=0; i<nrows_validation; i++)
    {
      S[ new_trainer_rows + i ] = Sx[ nrows_training + i ];
      E[ new_trainer_rows + i ] = Ex[ nrows_training + i ];
    }

  // complete indiv-IDs using above logic
  // i.e. we have not reset iid_idx, etc
  
  iid_start = -1;
  
  for (int i=nrows_training; i<nrows_training+nrows_validation; i++)
    {

      // only first row for this new person
      if ( iid_start == -1 ) iid_start = new_idx;	    
      iid_end = new_idx;
      ++new_idx;
    
      // advance to IID if at last epoch of current indiv?
      if ( i == Iendx[ iid_idx ] )
	{
	  // anything to add?
	  if ( iid_start != -1 )
	    {
	      Istart.push_back( iid_start );
	      Iend.push_back( iid_end );
	      I.push_back( Ix[iid_idx] );
	      // reset
	      iid_start = -1;
	    }
	  
	  // advance to next person 
	  if ( iid_idx < nix - 1 ) 
            ++iid_idx;
	}
    }

  
  // finally, update row count
  nrows_training = new_trainer_rows;

  // std::cout << "I szs = " << I.size() << " " << Istart.size() << " "<< Iend.size() << "\n";
  // std::cout << "X1 sxz " << X1.rows() << " " << E.size() <<" " << S.size() << "\n";
  
  // for (int i=0; i<I.size(); i++)
  //   std::cout << I[i] << "\t" << Istart[i] << "\t" << Iend[i] << "\n";
  
  // all done!
  
}


bool pops_t::attach_indiv_weights( const std::string & wlabel , bool training_dataset )
{
  
  // this is called at level 2 , prior to training only
  // we look for @vars ivars attached store, to find indiv-level weights
  // from the column named 'w' 
  
  // number of people w/ weights

  int cnt = 0;
  
  //
  // index people by the row at which they start in the data matrix
  // (to be passed to lgbm_t )
  //

  std::vector<uint64_t> starts;
  std::map<uint64_t,float> wtable;

  for (int i=0; i<Istart.size(); i++)
    {
      const bool in_training = holdouts.find( I[i] ) == holdouts.end () ;
      
      if ( ( training_dataset && in_training )
	   || ( (!training_dataset) && (!in_training) ) )	
	{
	  // get block start (potentially adjusting for offset into validation
	  uint64_t start = in_training ? Istart[i] : Istart[i] - nrows_training ;
	  starts.push_back( start );

	  // get weight
	  double w = 1;	  
	  if ( cmd_t::pull_ivar( I[i] , wlabel , &w ) ) ++cnt;
	  else w = 1; // not needed, to be make explicit: if missing keep w=1
	  wtable[ start ] = w;
	}
    }

  logger << "  updating weights for " << cnt << " of " << starts.size() 
	 << " individuals, from " << wlabel << " for the " 
	 << ( training_dataset ? "training" : "validation" ) << " dataset\n";

  //
  // make call to update LGBM weights
  //
  
  if ( training_dataset )
    lgbm.add_block_weights( lgbm.training , &lgbm.training_weights, starts, wtable );
  else
    lgbm.add_block_weights( lgbm.validation , &lgbm.validation_weights, starts, wtable );
	 
 
  return true;
}



bool pops_t::dump_weights() 
{
  std::string f = Helper::expand( pops_opt_t::model_weights_file );

  std::ofstream O1( f.c_str() , std::ios::out );

  logger << "  dumping weights to " << f << "\n";
  
  O1 << "ID\tTV\tSS\tW\n";
  
  std::vector<double> w1 = lgbm_t::weights( lgbm.training );
  std::vector<double> w2;
  if ( nrows_validation != 0 ) w2 = lgbm_t::weights( lgbm.validation );
  
  // std::cout << "w1, w2 size = " << w1.size() << " " << w2.size() << "\n";
  // std::cout << "X1 rows = " << X1.rows() << "\n";
  // std::cout << "nrow t, v = " << nrows_training << " " << nrows_validation << "\n";
  
  const int ni = I.size();
  int iid_idx = 0;

  //  std::cout << "I Iend = " << I.size() <<" " << Iend.size() << "\n";
  
  for (int i=0; i<X1.rows(); i++)
    {
      // training or validation?
      bool is_training = i < nrows_training ;
      
      // IID
      O1 << I[ iid_idx ] << "\t";
      
      O1 << ( is_training ? "T" : "V" ) << "\t";
      
      // epoch label (i.e. stage)
      O1 << pops_t::label( (pops_stage_t)S[i] ) << "\t";

      O1 << ( is_training ? i : i-nrows_training ) << "\t";

      O1 << ( is_training ? w1.size() : w2.size() ) << "\t";

      // weights
      O1 << ( is_training ? w1[i] : w2[i-nrows_training] ) << "\n";

      // advance IID if at last epoch of current indiv?
      if ( i == Iend[ iid_idx ] )
	{
	  if ( iid_idx < ni - 1 ) 
	    ++iid_idx;
	}
      
    }
  
  // all done
  O1.close();
  
  return true;
}


// helper to determine the presence of NaN in a matrix
// and construct a small report if so

pops_nan_report_t::pops_nan_report_t( const Eigen::MatrixXd & m )
{
  rows.clear(); 
  cols.clear();

  const bool anynans = m.array().isNaN().count() > 0 ;
  if ( ! anynans ) return;
  
  for (int c=0; c<m.cols(); c++) 
    {
      int cnt = m.col(c).array().isNaN().count();
      if ( cnt ) cols[c] = cnt;
    }

  for (int r=0; r<m.rows(); r++) 
    {
      int cnt = m.row(r).array().isNaN().count();
      if ( cnt ) rows[r] = cnt;
    }
  
}

bool pops_nan_report_t::any() const
{
  return rows.size() != 0 || cols.size() != 0;
}

void pops_t::stage_association1( Eigen::VectorXd & ftr , 
				 const std::vector<std::string> & ss )
{
  const int nobs = ftr.size();
  if ( nobs < 2 ) return;

  const int num_nan = ftr.array().isNaN().count();      
  writer.value( "N" , nobs - num_nan );
  writer.value( "N0" , num_nan );

  // no missing data 
  if ( num_nan == 0 ) 
    {         
      double pv = Statistics::anova( ss  , eigen_ops::copy_vector( ftr ) );
      if ( pv > -0.01 ) 
	{
	  if ( pv < 1e-200 ) pv = 1e-200;
	  writer.value( "ANOVA", -log10(pv) );
	}
      
      // nb. we have to standardize ftr first
      eigen_ops::scale( ftr , true , true );
      double wb = eigen_ops::between_within_group_variance( ss , ftr );
      writer.value( "WMAX", wb );
    }
  else
    {
      // need to splice out non-missing values
      const int nobs2 = nobs - num_nan ;
      if ( nobs2 < 2 ) return;
      Eigen::VectorXd ftr2 = Eigen::VectorXd::Zero( nobs2 );
      std::vector<std::string> ss2( nobs2 );
      int p = 0;
      for (int j=0; j<nobs; j++) 
	{
	  if ( Helper::realnum( ftr[j] ) )
	    {
	      ftr2[p] = ftr[j];
	      ss2[p] = ss[j];
	      ++p;
	    } 
	}
      
      double pv = Statistics::anova( ss2  , eigen_ops::copy_vector( ftr2 ) );
      if ( pv > -0.01 ) 
	{
	  if ( pv < 1e-200 ) pv = 1e-200;
	  writer.value( "ANOVA", -log10(pv) );
	}
      
      // nb. we have to standardize ftr first
      eigen_ops::scale( ftr2 , true , true );
      double wb = eigen_ops::between_within_group_variance( ss2 , ftr2 );
      writer.value( "WMAX", wb );

    }
  
}

void pops_t::stage_association() 
{
  // generate 1) a person x feature table
  //          2) an overall feature-level table

  logger << "  running feature/stage trainer associations\n";

  const int nobs = X1.rows();
  const int ncol = X1.cols();
    
  // feature labels
  std::vector<std::string> ftrs = pops_t::specs.select_labels();
  if ( ftrs.size() != ncol ) 
    Helper::halt( "interal error in stage_association()" );
  
  // stages
  std::vector<std::string> ss( nobs );
  for (int i=0; i<nobs; i++) 
    ss[i] = pops_t::label( (pops_stage_t)S[i] );
  
  //
  // consider each feature
  //
  
  for (int f=0; f<ftrs.size(); f++)
    {      
      writer.level( pops_t::specs.final2orig[ f ] + 1  , globals::feature_strat );
      Eigen::VectorXd ftr = X1.col(f);
      //std::cout << " OVERALL ftr " <<  f << "/" << ftrs.size() << "\n";
      stage_association1( ftr , ss );
    }
  writer.unlevel( globals::feature_strat );
  
  
  //
  // Now repeat for each individual
  //
    
  const int nt = Istart.size();

  for (int i=0; i<nt; i++)
    {
      std::cout << "i= " << i << " / " << nt << "\n";
      
      writer.level( I[i] , "TRAINER" );
		    
      // pull out data for this trainer only
      int fromi = Istart[i];
      int sz    = Iend[i] - Istart[i] + 1;

      std::cout << " fromi sz = " << fromi << "\t" << sz << "\n";

      Eigen::MatrixXd XI = X1.block( fromi , 0 , sz, ncol );
      
      // splice stages
      const int n1 = XI.rows();
      std::vector<std::string> ss1( n1 );
      for (int j=0; j<n1; j++) ss1[j] = ss[ fromi + j ];
      
      // consider each feature
      for (int f=0; f<ftrs.size(); f++)
	{	  
	  writer.level( pops_t::specs.final2orig[ f ] + 1 , globals::feature_strat );
	  Eigen::VectorXd ftr = XI.col(f);
	  //std::cout << " INDIV " << i << " ftr " <<  f << "/" << ftrs.size() << "\n";
	  stage_association1( ftr , ss1 );	  
	}
      writer.unlevel( globals::feature_strat );
    }
  writer.unlevel( "TRAINER" );  
  
}


#endif










