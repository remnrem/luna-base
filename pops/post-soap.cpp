
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
#include "stats/eigen_ops.h"
#include "stats/lda.h"
#include "miscmath/miscmath.h"

#include "db/db.h"
#include "helper/logger.h"

extern logger_t logger;
extern writer_t writer;


Eigen::MatrixXd pops_indiv_t::soap_X( bool * okay )
{

  //
  // Project full feature matrix (X1f) into 'nc' independent components for LDA (->U)
  //   i..e X1f == X1, but before any NaNs are set
  //
  
  const int nc = pops_opt_t::soap_nc ; 
  
  const int ne_full = P.rows();
    
  // find any columns w/ NaNs (e.g. un-specified covariates) and remove them
  std::set<int> included_Xcols;
  for (int c=0; c<X1f.cols(); c++)
    {
      bool okay = true;
      for (int r=0; r<X1f.rows(); r++)
	{
	  if ( std::isnan( X1f(r,c) ) )
	    {
	      okay = false;
	      break;
	    }	  
	}
      if ( okay ) included_Xcols.insert(c);
    }

  if ( included_Xcols.size() < 2 )
    {
      logger << "  ** could not find any X1 columns with non-missing values.. bailing on soap\n";
      *okay = false;
      Eigen::MatrixXd U;
      return U;
    }
  
  Eigen::MatrixXd X1ff = Eigen::MatrixXd::Zero( X1f.rows() , included_Xcols.size() );
  auto xx = included_Xcols.begin();
  int c = 0;
  while ( xx != included_Xcols.end() )
    {
      X1ff.col(c++) = X1f.col(*xx);
      ++xx;
    }
  
  Eigen::BDCSVD<Eigen::MatrixXd> svd( X1ff , Eigen::ComputeThinU | Eigen::ComputeThinV );
  Eigen::MatrixXd U = svd.matrixU().leftCols( nc );
  *okay = true;  
  return U;

}


double pops_indiv_t::simple_soap( const Eigen::MatrixXd & U ,
				  const std::vector<int> & ST )
{

  const int ne = ST.size();
  if ( U.rows() != ne ) Helper::halt( "internal error in pops_indiv_t::simple_soap()" );
  
  std::vector<std::string> sstr( ne );
  
  for (int e=0; e<ne; e++)
    sstr[e] = pops_t::labels5[ ST[e] ] ;

  lda_t lda( sstr  , U );
  
  lda_model_t lda_model = lda.fit();
    
  if ( ! lda_model.valid )
    {
      logger << "  *** could not fit SOAP model\n";
      return -1;
    }
  
  lda_posteriors_t prediction = lda_t::predict( lda_model , U ) ; 

  // get kappa
  double kappa = MiscMath::kappa( prediction.cl , sstr );

  return kappa;
  
}


void pops_indiv_t::apply_soap()
{
  
  //
  // Parameters
  //

  // SOAP-update confidence threshold
  double th = pops_opt_t::soap_threshold;
  
  // Stage-specific parameters
  // +5 : if stage has << 5 high confidence epochs, then leave as is
  // -5 : if stage has << 5 high confidence epochs, set all to missing (i.e. drop that stage)
  int mine = 5; 
  
  // stages with < mine confident epochs do not feature in the LDA
  //  when they decide:: are those re-assigned to the remaining confident epochs?
  //                or:: we leave them 'as is' 
  //   (could also only resassign if new conf > old conf or some such)
  
  const bool leave_rare_asis = true ; 
  
  //
  // number of PCs (for X1)
  //
  
  const int nc = pops_opt_t::soap_nc ; 


  //
  // Inputs : X1 -- data matrix
  //          P  -- posteriors
  //
  
  const int ne_full = P.rows();

  
  //
  // Project full feature matrix (X1f) into 'nc' independent components for LDA (->U)
  //   i..e X1f == X1, but before any NaN's set
  //


  // find any columns w/ NaN's (e.g. un-specified covariates) and remove them
  std::set<int> included_Xcols;
  for (int c=0; c<X1f.cols(); c++)
    {
      bool okay = true;
      for (int r=0; r<X1f.rows(); r++)
	{
	  if ( std::isnan( X1f(r,c) ) )
	    {
	      okay = false;
	      break;
	    }	  
	}
      if ( okay ) included_Xcols.insert(c);
    }

  if ( included_Xcols.size() < 2 )
    {
      logger << "  ** could not find any X1 columns with non-missing values.. bailing on soap\n";
      return;
    }
  
  Eigen::MatrixXd X1ff = Eigen::MatrixXd::Zero( X1f.rows() , included_Xcols.size() );
  auto xx = included_Xcols.begin();
  int c = 0;
  while ( xx != included_Xcols.end() )
    {
      X1ff.col(c++) = X1f.col(*xx);
      ++xx;
    }
  
  Eigen::BDCSVD<Eigen::MatrixXd> svd( X1ff , Eigen::ComputeThinU | Eigen::ComputeThinV );
  
  Eigen::MatrixXd U = svd.matrixU().leftCols( nc );

  logger << "  reducing " << X1f.cols() << " feature columns w/out NaN's to " << nc << " components\n";
    
  
  //
  // construct confidence and most-likely calls (w/ string labels)
  //  - already done, but just check here
  //
  
  std::vector<std::string> pops_predictions_str( ne_full );
  std::vector<int> pops_predictions( ne_full );
  std::vector<double> confidence( ne_full );
  
  for (int e=0; e<ne_full; e++)
    {
      int predx;
      double pmax = P.row(e).maxCoeff(&predx);
      // check prior work
      if ( predx != PS[e] ) Helper::halt( "book keeping error in POPS(1)" );
      confidence[e] = pmax;
      pops_predictions_str[e] = pops_t::labels5[ predx ] ;
      pops_predictions[e] = predx ;
    }
  
  
  
  //
  // Flag low-confidence assignments
  //

  std::vector<int> rows;
  std::vector<int> stg_count( P.cols() , 0 ); 
  
  for (int e=0; e< ne_full; e++)
    if ( confidence[e] >= th ) 
      {
	rows.push_back(e);
	stg_count[ pops_predictions[e] ]++;
      }
  
  int ne_conf = rows.size();
  
  
  //
  // Flag stages/classes that do not have enough unambiguous epochs
  //

  const int nstages_all = stg_count.size();
  
  std::set<int> low_conf_stages;
  std::vector<bool> stage_included( nstages_all );
  
  for (int s=0; s<nstages_all; s++)
    {
      stage_included[ s ] = stg_count[ s ] < mine ; 

      if ( stg_count[ s ] < mine )
	low_conf_stages.insert( s );
    }
  
  const int nstages_sufficient = nstages_all - low_conf_stages.size() ;

  if ( nstages_sufficient < 2 || ne_conf < 10 )
    {
      logger << "  ** fewer than two stages with a sufficient ( N > " << mine	     
	     << " ) number of unambiguous ( P > " << th << " ) epochs\n";
      logger << "  ** or less than 10 epochs with a confident call\n";
      return;
    }
  

  //
  // Identified rows to be spliced out
  //   (either low confidence, or could be high confidence, but for a stage
  //    with too few high confidence assignments)
  //
  
  ne_conf = 0;

  std::vector<bool> row_included( ne_full , true );
    
  for (int e=0; e<ne_full; e++)
    {
      if ( confidence[e] < th ||
	   low_conf_stages.find( pops_predictions[e] ) != low_conf_stages.end() )
	row_included[e] = false;
      else
	++ne_conf;
    }

  //
  // Make high-confidence sets for LDA
  //
  
  std::vector<std::string> S_conf( ne_conf );
  
  Eigen::MatrixXd U_conf = Eigen::MatrixXd::Zero( ne_conf , U.cols() );
  
  int r = 0;
  for (int e=0; e<ne_full; e++)
    {
      if ( row_included[e] )
	{
	  S_conf[ r ] = pops_predictions_str[ e ];
	  U_conf.row( r ) = U.row( e );
	  ++r;
	}      
    }


  // std::cout << "S_conf \n";
  
  // for (int l=0; l<S_conf.size(); l++)
  //   {
  //     std::cout << S_conf[l] ;
  //     for (int q=0; q<U_conf.cols(); q++)
  // 	std::cout << "\t" << U_conf(l,q) ;
  //     std::cout << "\n";
  //   }

  //
  // LDA-based SOAP on unambiguous epochs 
  //
  
  lda_t lda( S_conf  , U_conf );     

  // get 'S' mapping
  
  //  std::
  
  
  // note: second 'S" param means to set priors based on the full/original st[] rather
  // than the subset of unambiguous values
  
  //lda_model_t lda_model = lda.fit( false , &pops_predictions_str );
  lda_model_t lda_model = lda.fit( false );
  
  // std::cout << "grops means = \n" << lda_model.means << "\n\n";
  // std::cout << "grops scaling = \n" << lda_model.scaling << "\n\n";
  
  if ( ! lda_model.valid )
    {
      logger << "  *** could not fit SOAP model, leaving posteriors unaltered\n";
      return;
    }
  
  
  //
  // Predict back all rows using this model
  //
  
  lda_posteriors_t prediction = lda_t::predict( lda_model , U ) ; 
  
  // labels in lda_model.labels[]
  std::vector<int> old2new( 5, -1 );
  for (int i=0; i<lda_model.labels.size(); i++)
    {
      if      ( lda_model.labels[i] == "W" ) old2new[0] = i;
      else if ( lda_model.labels[i] == "R" ) old2new[1] = i;
      else if ( lda_model.labels[i] == "N1" ) old2new[2] = i;
      else if ( lda_model.labels[i] == "N2" ) old2new[3] = i;
      else if ( lda_model.labels[i] == "N3" ) old2new[4] = i;      
    }
  

    //
  // Dumper
  //

  if ( 0 )
    {

      for (int e=0; e < ne_full; e++)
	{
	  
	  std::cout << e << "\t";
	  
	  // original prediction
	  std::cout << pops_predictions_str [e] << "\t" ;

	  // was this included?
	  std::cout << row_included[ e ] << "\t";

	  // original PP (W,R,N1,N2,N3)
	  for (int j=0; j < nstages_all ; j++)
	    std::cout << " " << P(e,j);
	  
	  std::cout << " -->\t ";

	  // if high-confidence, keep as is
	  if ( row_included[ e ] )
	    {
	      
	      std::cout << ".\t.";
	      
	      // for (int j=0; j < nstages_all ; j++)
	      // 	std::cout << "\t" << P(e,j);	      
	      
	    }
	  else
	    {
	      // this is an ambiguous epoch --> what does the new model say?
	      
	      std::cout << "\t" << prediction.cl[ e ]
			<< "\t" << ( pops_predictions_str[ e ]  != prediction.cl[ e ] ? "X" : "." );

	      // use 'old2new[] to get order; -1 if not included
	      for (int j=0; j<5; j++)
		{
		  if ( old2new[j] == -1 )
		    std::cout << " .";
		  else
		    std::cout << " " << prediction.pp( e , old2new[j] );
		  
		}
	    }
	  
	  std::cout << "\n";
	}
    }
  
  
  //
  // Modify originals
  //
  
  int nchanged = 0;
  
  for (int e=0; e < ne_full; e++)
    {
      
      const bool low_conf = ! row_included[ e ] ; 
      
      const bool retain = leave_rare_asis && low_conf_stages.find( PS[e] ) != low_conf_stages.end() ;
      
      if ( low_conf && ! retain )
	{

	  // original prediction:: PS[e]
	  // original posteriors:: P(e,j);
	  //                       confidence[e]
	  
	  // revised prediction:: prediction.cl[ e ]  [ str ]
	  
	  int revised = 0; // W
	  if      ( prediction.cl[ e ] == "R" ) revised = 1;
	  else if ( prediction.cl[ e ] == "N1" ) revised = 2;
	  else if ( prediction.cl[ e ] == "N2" ) revised = 3;
	  else if ( prediction.cl[ e ] == "N3" ) revised = 4;
	  
	  double revised_conf = P.row(e).maxCoeff();

	  // only update if the new conf is greater than the old one
	  // nb. we may have fewer classes here, so not exactly apples-to-apples
	  
	  if ( revised_conf > confidence[e] )
	    {
	      
	      for (int j=0; j<5; j++)
                {
                  if ( old2new[j] == -1 )
                    P(e,j) = 0;
                  else
                    P(e,j) = prediction.pp( e , old2new[ j ] );
		}
	      
	      // update most likely predicted stage
	      PS[e] = revised;
	      if ( PS[e] != revised ) ++nchanged;
	    }
	}
    }
  
  //
  // all done
  //

  logger << "  changed " << nchanged << " epochs based on soap cleaning\n";
    
  
}



void pops_indiv_t::grid_soap()
{

  //
  // iteratively circle through stages, doing grid search to optimize R
  //

  
  // likelihood rescaling factors:: initially all at 1.0, i.e. no changes
  
  Eigen::VectorXd r = Eigen::VectorXd::Constant( 5 , 1 );
  
  // make SOAP feature matrix

  bool okay = true;
  Eigen::MatrixXd U = soap_X( &okay );
  if ( ! okay ) return;
  
  // baseline SOAP 

  double k = simple_soap( U , PS );

  // update PS given P, count # of stages
  std::vector<int> cnts0(5,0);  
  int nstages_orig = update_predicted(&cnts0);  
  
  //
  // grid search around rescaled likelhoods, and use SOAP kappa to determine the optimal values
  //
  
  // track max kappa 
  double max_kappa = 0;
 
  // if global priors have not been specified (i.e. if not running es-priors) then
  // set to uniform - i.e. should not matter as we are not changing them

  if ( pops_t::ES_global_priors.size() == 0 )
    pops_t::ES_global_priors = Eigen::VectorXd::Constant( 5, 0.2 );
  
  // copy posteriors

  Eigen::MatrixXd P0 = P;

  std::vector<double> ll = MiscMath::linspace( pops_opt_t::lk_lwr , pops_opt_t::lk_upr , pops_opt_t::lk_steps );


  //
  // BY default, only rescale likelihoods for stages that are ~never otherwise confidently 
  // assigned.  By this on % of epochs that have a confidence above X 
  //
  
  std::vector<int> stgs;
  const double conf_th1 = pops_opt_t::soap_grid_mean_conf; // 0.8 by default
  const double conf_th2 = 0.05;
  // i.e. if more than 5% of epochs have P(stage|data) > 0.8, then do not try to 
  // rescale this stage, as it already has sufficient confident assignments
  
  for (int ss=0; ss<5; ss++)
    {
      double prop = 0 ; 
      const int ne = P.rows();
      for (int e=0; e<ne; e++)
	if ( P(e,ss) > conf_th1 ) ++prop; 
      prop /= (double)ne;
      if ( prop < conf_th2 )
	{
	  stgs.push_back( ss );
	  logger << "  SOAP-scaling likelihoods for " 
		 << pops_t::labels5[ ss ] << " ( " << prop << " epochs > " << conf_th1 << " conf )\n";
	}
    }
  
  //
  // rescale stages identified above
  //
  
  for (int ss=0; ss<stgs.size(); ss++)
    {
      int s2 = stgs[ss] ; 
      
      //logger << "  rescaling " << pops_t::labels5[ s2 ] << " likelihoods\n";
      writer.level( pops_t::labels5[ stgs[ss] ] , globals::stage_strat ); 
            
      double max_fac = 1;
  
      for (double sidx = 0; sidx < ll.size(); sidx++)
	{
	  
	  //
	  // likelihood rescaling factors
	  //
	  
	  double fac = ll[sidx];
	  
	  //Eigen::VectorXd r2 = Eigen::VectorXd::Constant( 5 , 1 );
	  //r2[ss] = fac;

	  // alter this stage
	  r[ s2 ] = fac;
	  
	  // set to the originals
	  P = P0;

	  // rescale posteriors
	  const int ne = P.rows();
	  
	  for (int e=0; e<ne; e++)
	    P.row(e) = update_posteriors( P.row(e) ,
					  pops_t::ES_global_priors ,
					  NULL , // no change in priors
					  &r );
	  
	  // update PS given P
	  std::vector<int> cnts1(5,0);
	  int nstages = update_predicted(&cnts1);
	  
	  // redo SOAP (w/ same U)	  
	  double k1 = nstages >= nstages_orig ? simple_soap( U , PS ) : 0 ; 
	  
	  writer.level( fac , "FAC" );
	  writer.value( "K" , k1 );
	  writer.value( "NS" , nstages );

	  if ( 0 )
	    {
	      std::cout << s2 << "\t" << fac << "\t"
			<< nstages << "/" << nstages_orig << "\t"
			<< max_kappa << "\t" << k1
			<< " r= " << r.transpose() << " :";
	      
	      for (int i=0; i<5; i++) std::cout << " " << cnts0[i];
	      std::cout <<" || ";
	      
	      for (int i=0; i<5; i++) std::cout << " " << cnts1[i];
	      std::cout <<"\n";
	    }
	  
	  
	  if ( k1 > max_kappa )
	    {
	      //std::cout << " UPDATING k = " << k1 << "\n";
	      max_kappa = k1;
	      max_fac = fac;
	    }
	  
	}

      writer.unlevel( "FAC" );

      
      //
      // Final rescaling for this stage
      //
      
      P = P0;
      
      r[ s2 ] = max_fac;
      
      for (int e=0; e<P.rows(); e++)
	P.row(e) = update_posteriors( P.row(e) ,
				      pops_t::ES_global_priors ,
				      NULL , 
				      &r );
      update_predicted();
      
      //logger << "  reacaled likelihood R = " << r.transpose() << "\n";

      writer.value( "RESCALE_REM_FAC" , max_fac );
      writer.value( "RESCALE_REM_K0" , k );
      writer.value( "RESCALE_REM_K1" , max_kappa );
      
      // next stage (curr. REM only)
    }

  writer.unlevel( globals::stage_strat );
}




#endif
