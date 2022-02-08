
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

#include "stats/lgbm.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "stats/eigen_ops.h"
#include "db/db.h"
#include "dsp/tv.h"

extern logger_t logger;
extern writer_t writer;

pops_opt_t pops_t::opt;
lgbm_t pops_t::lgbm;
pops_specs_t pops_t::specs;

  
//
// create a level 2 feature library and save
// i.e. for trainer and/or validation library
//

void pops_t::make_level2_library( param_t & param )
{

  std::string lgbm_model  = param.requires( "model" );
  std::string data_file   = param.requires( "data" );

  std::string lgbm_config   = param.has( "config" ) ? param.value( "config" ) : ".";
  std::string feature_file  = param.has( "features" ) ? param.value( "features" ) : "." ;
  
  // feature specifications (used to generate l1-data)
  pops_t::specs.read( feature_file );

  // get previous data: assume single file, concatenated
  // this will populate: X1, S, E and Istart/Iend
  load1( data_file );

  // expand X1 to include space for level-2 features
  X1.conservativeResize( Eigen::NoChange , pops_t::specs.na );
  
  // std::cout << " Is " << Istart.size();
  // for (int i=0; i<Istart.size(); i++)
  //   std::cout << "i " << i+1 << " " << Istart[i] << " ... " << Iend[i] << "\n";

  const int n = S.size();
  std::map<int,int> ss;
  for (int i=0; i<n; i++) ss[S[i]]++;
  std::map<int,int>::const_iterator ii = ss.begin();
  while ( ii != ss.end() )
    {
      logger << "  " << pops_t::label( (pops_stage_t)ii->first ) << " = " << ii->second << "\n";
      ++ii;
    }
  
  // derive level-2 stats
  level2();

  // LGBM config (user-specified, or default for POPS)
  if ( lgbm_config == "." )
    lgbm.load_pops_default_config();
  else    
    lgbm.load_config( lgbm_config );

  // do training
  fit_model( lgbm_model );
  
}


//
// derive level 2 stats (from pops_t::specs)
//

void pops_t::level2()
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
      
      logger << "  level-2 features: " << l2ftr << " ( " << from_block << " --> " << to_block << " )\n";

      logger << "  from cols (N=" << from_cols.size() << ")\n";
      for (int j=0;j<from_cols.size();j++) logger << " " << from_cols[j] ;
      logger << "\n";
      
      logger << "  to cols (N=" << to_cols.size() << ")\n";
      for (int j=0;j<to_cols.size();j++) logger << " " << to_cols[j] ;
      logger << "\n\n";

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
      // SVD
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
	  
	  // SVD
	  Eigen::BDCSVD<Eigen::MatrixXd> svd( D , Eigen::ComputeThinU | Eigen::ComputeThinV );
	  Eigen::MatrixXd U = svd.matrixU();
	  Eigen::MatrixXd V = svd.matrixV();
	  Eigen::VectorXd W = svd.singularValues();

	  // copy U back to X1
          for (int j=0; j<nto; j++) X1.col( to_cols[j] ) = U.col(j);
	  
	  // save W and V to a file (i.e. for use later in prediction models, to project
	  // test cases into this space)
	  std::ofstream OUT1( Helper::expand( wvfile ).c_str() , std::ios::out );
	  OUT1 << V.rows() << " " << nc << "\n";
	  for (int i=0;i<V.rows(); i++)
	    for (int j=0;j<nc; j++)
	      OUT1 << " " << V(i,j) ;
	  OUT1 << "\n";
	  for (int j=0;j<nc; j++)
	    OUT1 << " " << W[j] ;
	  OUT1 << "\n";
	  OUT1.close();
	  
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

  if ( 1 )
    {
      std::ofstream OUT2( "X1.mat" , std::ios::out );
      OUT2 << X1 << "\n";
      OUT2.close();
    }
  
}


//
// fit and save a LGBM model
//

void pops_t::fit_model( const std::string & modelfile )
{
  std::cout << "S1\n";
  // trainging data
  lgbm.attach_training_matrix( X1 );

  std::cout << "S2\n";
  // labels
  lgbm.attach_training_labels( S );

  std::cout << "S3\n";

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

  lgbm.load_config( param.requires( "config" ) );
  lgbm.load_model( param.requires( "model" ) );

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

#endif
