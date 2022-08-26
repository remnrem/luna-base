
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

#include "db/db.h"
#include "helper/logger.h"

extern logger_t logger;
extern writer_t writer;

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



#endif
