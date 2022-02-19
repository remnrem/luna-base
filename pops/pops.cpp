
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
#include "helper/helper.h"
#include "helper/logger.h"
#include "stats/eigen_ops.h"
#include "db/db.h"
#include "dsp/tv.h"

#include "helper/zfile.h"

extern logger_t logger;
extern writer_t writer;

//pops_opt_t pops_t::opt;
lgbm_t pops_t::lgbm;
bool pops_t::lgbm_model_loaded = false;
pops_specs_t pops_t::specs;

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

  // always requires data input
  const std::string data_file   = param.requires( "data" );

  // needs a feature file, but will use default/internal ('.') if one is not specified
  const std::string feature_file  = param.has( "features" ) ? param.value( "features" ) : "." ;

  //
  // *either*, dump features after constructing feature matrix, or run LGBM
  //

  std::string lgbm_model = "";
  std::string lgbm_config = "";
    
  const bool dump_feature_matrix = param.has( "dump" );
  bool is_trainer = true;
  std::string dump_file = "";

  if ( dump_feature_matrix ) 
    {      
      if ( param.has( "model" ) || param.has( "config" ) )
	Helper::halt( "cannot specify both dump and model/config" );
      if ( param.value( "dump" ) == "test" )
	is_trainer = false;
      else if ( param.value("dump" ) != "training" )
	Helper::halt( "'dump' must be set to 'training' or 'test'" );
      dump_file = param.requires( "file" );
    }
  else 
    {
      // else we need LGBM configs - always needs a model
      lgbm_model = param.requires( "model" );
      // can have a default config file
      lgbm_config = param.has( "config" ) ? param.value( "config" ) : ".";
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

  //  
  // get previous data: assume single file, concatenated
  // this will populate: X1, S, E and Istart/Iend
  // validation IDs will be at the end
  //

  load1( data_file );

  //
  // expand X1 to include space for level-2 features
  //

  X1.conservativeResize( Eigen::NoChange , pops_t::specs.na );

  
  // 
  // summarize training and validation dataset counts  
  //

  std::map<int,int> ss_training, ss_valid;
  for (int i=0; i<nrows_training; i++) 
    ss_training[S[i]]++;
  for (int i=nrows_training; i<nrows_validation+nrows_training; i++) 
    ss_valid[S[i]]++;
  
  logger << "  in nT=" << Istart.size() - holdouts.size() << " training individuals, and nV=" 
	 << holdouts.size() << " validation individuals, stage epoch counts:\n";
  std::map<int,int>::const_iterator ii = ss_training.begin();
  while ( ii != ss_training.end() )
    {
      logger << "  " << pops_t::label( (pops_stage_t)ii->first ) << "\t" 
	     << " train = " << ii->second << "\t"
	     << " validation = " << ss_valid[ ii->first ] << "\n";
      ++ii;
    }

  //  
  // derive level-2 stats (jointly in *all* indiv. training + valid ) 
  //  -- typically trainer, but allow for test individuals here if we
  //     are using the --pops command not to train, but with dump=test
  //     (this will instruct level2() to read the SVD rather than create)
  //
  
  level2( is_trainer );
  
  
  //
  // output the the feature matrix 
  //

  if ( dump_feature_matrix ) 
    {
      dump_matrix( dump_file );
      return;
    }


  //
  // LGBM config (user-specified, or default for POPS)
  //

  if ( lgbm_config == "." )
    lgbm.load_pops_default_config();
  else    
    lgbm.load_config( lgbm_config );

  // set number of iterations?
  if ( param.has( "iterations" ) ) 
    lgbm.n_iterations = param.requires_int( "iterations" );

  //  
  // stage weights?    
  //

  std::vector<double> wgts(  pops_opt_t::n_stages , 1.0 );
  
  if ( param.has( "weights" ) ) // order must be W, R, NR...
    {
      wgts = param.dblvector( "weights" );
      if ( wgts.size() != pops_opt_t::n_stages )
	Helper::halt( "expecting " + Helper::int2str( pops_opt_t::n_stages ) + " stage weights" );
    }
  
  lgbm_label_t weights( pops_opt_t::n_stages == 5 ? pops_t::labels5 : pops_t::labels3 , wgts );

  //
  // do training
  //
 
  fit_model( lgbm_model , weights );
  
}


//
// derive level 2 stats (from pops_t::specs)
//

void pops_t::level2( const bool training )
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
      
      logger << "   - adding level-2 feature " << l2ftr << ": "; 

      if ( from_cols.size() > 0 )  // time track does not have "from" cols
	logger << from_block << " (n=" << from_cols.size() << ") ";
      
      logger << "--> " 
	     << to_block << " (n=" << to_cols.size() << ", cols:" 
	     << to_cols[0] << "-" << to_cols[to_cols.size()-1] << ") \n";
      
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
		  D.segment( fromi , sz ) = eigen_ops::moving_average( D.segment( fromi , sz ) , fwin );
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
      // SVD  - done differently for trainers versus targets
      //
      
      if ( spec.ftr == POPS_SVD )
	{

	  const int nc = spec.narg( "nc" );
	  const std::string wvfile = spec.arg[ "file" ];
	  
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
	      
	      logger << "   - writing SVD W and V to " << wvfile << "\n";

	      std::ofstream OUT1( Helper::expand( wvfile ).c_str() , std::ios::out );
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
		  // allow for 'path' option to modify where this file is
		  std::string filename = pops_t::update_filepath( wvfile );
		  logger << "   - reading SVD W and V from " << filename << "\n";
		  if ( ! Helper::fileExists( filename ) ) 
		    Helper::halt( "cannot find " + filename + "\n (hint: add a 'path' arg to point to the .svd file" );
		  std::ifstream IN1( filename.c_str() , std::ios::in );
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

void pops_t::fit_model( const std::string & modelfile , const lgbm_label_t & weights )
{
  // training   X1.topRows( nrows_training )   -->   S1
  // validaiton X1.bottomRows( nrows_validation )   --> S2

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
  lgbm.apply_label_weights( lgbm.training , weights );
  
 
  // validation data?
  if ( nrows_validation ) 
    {
      lgbm.attach_validation_matrix( X1.bottomRows( nrows_validation) );

      // labels
      lgbm.attach_validation_labels( S2 );
      
      // training weights
      lgbm.apply_label_weights( lgbm.validation , weights );
      
    }

  // fit model
  lgbm.create_booster();
  
  // all done, write model
  lgbm.save_model( modelfile );
  
}



//
// attach a final LGBM model
//

void pops_t::load_model( param_t & param )
{
  if ( ! lgbm_model_loaded )
    {
      lgbm.load_config( param.requires( "config" ) );
      lgbm.load_model( param.requires( "model" ) );
      lgbm_model_loaded = true;  
    }
}

void pops_t::outliers( const Eigen::VectorXd & x ,
		       const double th ,
		       const std::vector<int> & staging , 
		       std::vector<int> * staging2 )
{
  double sum = 0 , sumsq = 0;
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
  const double mean = sum/n;
  const double meansq = mean * mean ; 
  const double sd = sqrt( sumsq / (n-1) - ( n/(n-1) ) * meansq );
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
  //   0 all epochs   A
  //   1 only epochs with similar flanking observed stages (i.e. 'consistent' sleep)  A-A-A
  //   2 only left-epochs at a transition (i.e. if the following obs epoch is not the same)  A-B
  //   3 only right-epochs at a transition (i.e. if the prior obs epoch was not the same)  B-A
  //   4 only 'singleton' epochs flanked by the same epoch on both sides  B-A-B
  //   5 only 'singleton' epochs, with any flanking epochs  B-A-C

  //  further, if ostage != -1, then oonly lookat epochs with that obs stage type
  
  std::vector<int> obs; 
  std::vector<int> pred;
  
  const int ne = obs_.size();

  if ( type == 0 ) 
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
	  else if ( type == 2 ) // *-A-B 
	    add = right_disc;
	  else if ( type == 3 ) // B-A-*
	    add = left_disc;
	  else if ( type == 4 ) // B-A-B
	    add = left_disc && ! left_right_disc ;	    
	  else if ( type == 5 ) // B-A-C
	    add = left_disc && left_right_disc ;
	  
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

      logger << "\t   Obs:";
      std::set<int>::const_iterator uu = uniq.begin();
      while ( uu != uniq.end() )
	{
	  logger << "\t" << pops_t::label( (pops_stage_t)(*uu) );
	  ++uu;
	}
      logger << "\tTot\n";	
      
      logger << "  Pred:";
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

      
      // conditional probabilties  / res[][] / cols[] 
      uu = uniq.begin();
      while ( uu != uniq.end() )
        {
	  writer.level( pops_t::label( (pops_stage_t)(*uu) ) , "PRED" );
	  std::set<int>::const_iterator jj = uniq.begin();
          while ( jj != uniq.end() )
	    {
	      writer.level( pops_t::label( (pops_stage_t)(*jj) ) , "OBS" );
	      writer.value( "N" , res[ *uu ][ *jj ] );
	      if ( cols[ *uu ] > 0 ) 
		writer.value( "P" , res[ *uu ][ *jj ] / cols[ *jj ] );
	      ++jj;
	    }
	  writer.unlevel( "OBS" );
	  ++uu;
	}
      writer.unlevel( "PRED" );
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
  if ( f2[0] != globals::folder_delimiter && pops_opt_t::pops_path != "" )
    f2 = globals::folder_delimiter + pops_opt_t::pops_path + f2;
  
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
      Z1 << pops_t::label( (pops_stage_t)S[i] );
      for (int j=0; j< X1.cols(); j++)
	Z1 << "\t" << X1(i,j);
      Z1 << "\n";
    }
  Z1.close();
}


#endif






