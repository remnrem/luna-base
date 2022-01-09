
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


// RESOAP : iterative staging 

void suds_indiv_t::resoap_alter1( edf_t & edf , int epoch , suds_stage_t stage )
{
    
  // actual number of epochs
  
  const int ne_all = edf.timeline.num_epochs();
  
  // for SOAP, number of included epochs based on good signal data 
  // (i.e. as all epochs are included for SOAP targets, irrespective for 
  // ovserved stage being known or not)
  
  const int ne_included = obs_stage_valid.size();
  
  // nb. 'epoch' is 1-based , given by the user
  if ( epoch < 1 || epoch > ne_all ) 
    Helper::halt( "bad epoch value, outside range" );
  
  // some epochs may be skipped, e.g. due to signal outliers
  
  // valid epochs  : std::vector<int> epochs;   
  // all stages    : std::vector<suds_stage_t> obs_stage; 
  // valid stages  : std::vector<suds_stage_t> obs_stage_valid; 
  // same, but str : std::vector<std::string> y;   (send to lda_t() 
  
  // we need to update only y[i]; so the 'original' is kept in obs_stage(_valid)
    
  //
  // Update 'y' and check we have en
  //
  
  // epochs[] contains the codes of epochs actually present in the model/valid                                                                                                                      
  //  (in 0, ... encoding)
  
  bool updated = false;
  
  for (int i=0; i < epochs.size(); i++) 
    {
      
      // internal epoch number = epochs[i]
      // display epoch = edf.timeline.display_epoch( i )
      // nb. in SOAP context, with no restructuring of the EDF, the display epoch
      // will typically be +1 the internal epoch, i.e. not expecting discontinous
      // codes;  the user is expected to 
      
      // for y and obs_stage_valid
      int e0 = i;
      
      // for obs_stage
      int e1 = epochs[i];
      
      // for user disokay epoch (1-based)
      int e2 = edf.timeline.display_epoch( e1 );

      // update this single 'observed' stage
      if ( epoch == e2 ) 
	{
	  logger << "  changing epoch " << epoch << " from " << y[e0] << " to " << suds_t::str( stage ) << "\n";
	  y[e0] = suds_t::str( stage );
	  // obs_stage_valid[ e0 ] = stage;
	  // obs_stage[ e1 ] = stage;
	  updated = true;
	}
      
      // track what we have     
    }
  
  if ( ! updated ) 
    logger << "  no updates made: did not find epoch " << epoch << " (with valid signal data)\n";
  
  
}


void suds_indiv_t::resoap_pickN( edf_t & edf , int pick )
{
  // for evaluation of SOAP only: 

  // pick N each of 'labels', using the original observed stages
  if ( obs_stage_valid.size() != y.size() )
    Helper::halt( "cannot use RESOAP pick without original staging" );
  
  // first scrub
  for (int i=0; i < suds_t::cached.y.size(); i++)
    suds_t::cached.y[i] = suds_t::str( SUDS_UNKNOWN );
  
  const int nss = suds_t::labels.size();

  std::map<std::string,int> scounts;

  // N or more  versus exactly N
  bool exact = pick < 0 ;
  if ( exact ) pick = -pick;

  const int n = y.size();

  // Yates-Fisher shuffle to get a random ordering
  std::vector<int> a( n );
  CRandom::random_draw( a );
  
  std::set<std::string> done;
  for (int i=0; i<n; i++)
    {
      
      int p = a[i]; // random draw

      std::string ss = suds_t::str( obs_stage_valid[p] );
      if ( ss == "?" ) continue;
	       
      if ( exact )
	{
	  // only add if fewer than needed?
	  int c = scounts[ ss ];	  
	  if ( c < pick )
	    {
	      y[p] = ss;
	      ++scounts[ ss];	      
	    }
	}
      else
	{	    
	  y[p] = ss;
	  ++scounts[ ss ];
	}
      
      // done for this stage?
      if ( scounts[ y[p] ] == pick )
	done.insert( y[p] );

      // all done?
      if ( done.size() == nss ) break;	
    }

}


void suds_indiv_t::resoap( edf_t & edf , bool epoch_level_output )
{

  logger << "  re-SOAPing...\n";

  //
  // this impacts format of epoch-level output
  //

  suds_t::soap_mode = 2;
  
  //
  // Count "observed" stages
  //
  
  const int n = y.size();
  
  std::map<std::string,int> ycounts;
  for (int i=0; i<n; i++) ++ycounts[ y[i] ];
  

  //
  // requires at leasrt two stages w/ at least 3 observations, and has to 
  // be greater than the number of PSCs
  //
  
  const int required_n = 3;
  
  logger << "  epoch counts:";
  int s = 0;
  int t = 0;
  int tt = 0;
  std::map<std::string,int>::const_iterator yy = ycounts.begin();
  while ( yy != ycounts.end() )
    {
      logger << " " << yy->first << ":" << yy->second;     
      tt += yy->second;
      if ( yy->first != "?" && yy->second >= required_n ) 
	{
	  ++s;
	  t += yy->second;
	}
      ++yy;
    }
  logger << "\n";

  writer.value( "S" , s );
  writer.value( "OBS_N" , t ); // need at least 3 of each for 't'
  writer.value( "OBS_P" , t/(double)tt );

  
  bool okay = s >= 2 ;
  
  // for p predictors, require at least p+2 observations
  if ( ! ( t > nc+1 ) ) okay = false;
  
  if ( ! okay )
    {
      logger << "  not enough non-missing stages for LDA with " << nc << " predictors\n";
      writer.value( "FIT" , 0 );
      return;
    }
  
  
  //
  // Re-fit the LDA
  //
  
  Eigen::MatrixXd pp;
  
  int dummy = self_classify( NULL , &pp );

  if ( dummy == 0 )
    {
      logger << "  LDA model could not converge with the current stage proposal\n";
      writer.value( "FIT" , 0 );      
      return;
    }


  //
  // Model okay
  //

  writer.value( "FIT" , 1 );


  //
  // output stage probabilities 
  //
  
  const double epoch_sec = edf.timeline.epoch_length();
  
  std::vector<std::string> final_pred = suds_t::max( pp , lda_model.labels );
  
  summarize_kappa( final_pred , true );
  
  // actual number of epochs
  const int ne_all = edf.timeline.num_epochs();

  const int bad_epochs = summarize_stage_durations( pp , lda_model.labels , ne_all , epoch_sec );

  if ( epoch_level_output )
    summarize_epochs( pp , lda_model.labels , ne_all , edf );
  
}


// struct suds_pp_sorter_t {
//   suds_pp_sorter_t( double p , int e ) : p(p) , e(e) { } 
//   double p;
//   int e;
//   bool operator<( suds_pp_sorter_t & rhs ) const
//   {
//     // flip, to get highest first
//     if ( p > rhs.p ) return true;
//     if ( p < rhs.p ) return false;
//     return e < rhs.e;
//   }
// };

int suds_indiv_t::resoap_update_pp( std::vector<std::string> * st , 
				    Eigen::MatrixXd * pp ,
				    const std::vector<std::string> & labels , 
				    const bool global_mode ) 
{

  //
  // Inputs:: proposal staging       st[] - all epochs will be non-missing
  //          posterior probs        pp[,]
  //          labels for posteriors  labels
  //          global_mode            allow different params based on condition

  //          target.U
  
  //          suds_t::soap_[global]_update_th            soap=0.8,5  soap1=0.8
  //          suds_t::soap_[global]_update_min_epochs     or soap=0.9,-5 

  //
  // Outputs  -- modified stages (st) and posterior (pp) based on SOAP
  //          -- ret val = number of epochs changed
  
  //  SOAP model only inlcudes 'ambiguously' assigned epochs (i.e. above thresholds)
  //          -- for epochs without ambiguous assignment:::
  //            1) if most likely stage has enough (e.g. min_epochs = 5 above) unambiguous epochs, we just set to ? (so SOAP fills in)
  //            2a) otherwise, we *either* leave those as is, and only do SOAP on the other stages (minus those cols)
  //            2b) OR, we drop that entire class label, and do SOAP for all epochs  
  //            3) OR, if min_epochs == 0 
    
  // global_mode :  whether called for the final set (T)
  //                or per trainer individual prediction (F)
  //    i.e. we might want to treat edge cases differently (whether to drop things, etc)

  
  const double th   = global_mode ? suds_t::soap_global_update_th : suds_t::soap_update_th ;
  const int    mine = global_mode ? abs( suds_t::soap_global_update_min_epochs ) : abs( suds_t::soap_update_min_epochs );
  const bool   leave_rare_asis = global_mode ? suds_t::soap_global_update_min_epochs > 0 : suds_t::soap_update_min_epochs > 0;
  
  const int rows = pp->rows();
  const int cols = pp->cols();
  
  if ( st->size() != rows ) 
    Helper::halt( "internal error in resoap_update_pp(), w/ rows" );
  std::cout << " labels.siize + " << labels.size() << " " << cols << "\n";
  if ( labels.size() != cols )
    Helper::halt( "internal error in resoap_update_pp(), w/ cols" );

  // counts: (high-conf) stages, epochs 
  std::set<std::string> stgs;
  std::map<std::string,int> stgs2;

  // track most confident epochs within each class
  //std::map<std::string,std::set<suds_pp_sorter_t> > stgpp;

  // make a copy of st, which (below) may have epochs flagged as unknown if ambiguous
  std::vector<std::string> st2 = *st; 
    
  //
  // flag low-confidence assignments
  //
  
  int blanked = 0;

  for (int i=0; i<rows; i++)
    {
      
      const double mx = pp->row(i).maxCoeff();
      
      if ( mx >= th ) 
	stgs2[ st2[i] ]++;
      else // this epoch was not 'unambiguous'
	{
	  ++blanked;
	  st2[i] = suds_t::str( SUDS_UNKNOWN );	  
	}

      // also, total number of stages seen      
      stgs.insert( st2[i] );

      // also, track best pp from each stage
      // for (int j=0;j<cols;j++)	
      // 	stgpp[ st2[i] ].insert( suds_pp_sorter_t( (*pp)(i,j) , i ) ); 
      
    }
  
  
  //
  // Flag stages/classes that do not have enough unambiguous epochs
  //
  
  std::vector<bool> col_included( cols , true );
  std::set<std::string> asis;
  int drop_cols = 0;
  
  for (int s=0; s<cols; s++)
    if ( stgs2[ labels[s] ] < mine )
      {
	col_included[s] = false;
	asis.insert( labels[s] );
	drop_cols++;
      }

  //
  // Any columns flagged above will be dropped; the question here is
  // whether the rows containing those stages (as most-likely
  // assignment) will also be excluded from SOAP (i.e. left as-is, and
  // keep that stage in the final answer) or whether we will let them
  // be assigned to another most likely class (i.e. drop that stage
  // from the final answer)
  //
  
  std::vector<bool> row_included( rows , true );
  int drop_rows = 0;

  if ( leave_rare_asis )
    {      
      for (int i=0; i<rows; i++)
	{
	  // if we want to leave this class as is, remove the column and then
	  // put the original labels back into stg2; these will be spliced back
	  // into the final, adjusted staging
	  
	  if ( asis.find( (*st)[i] ) != asis.end() )
	    {
	      // copy back original 
	      st2[i] = (*st)[i] ;

	      // but flag to drop this from SOAP
	      row_included[i] = false;
	      ++drop_rows;
	    }
	}
    }
  
  //
  // Splice out columns and rows that will not go into SOAP
  //

  Eigen::MatrixXd U2 = U;
  
  std::vector<std::string> S = st2;
  
  //
  // Splice out any rare epochs (i.e. leave those as is)
  //
  
  if ( drop_cols || drop_rows )
    {
      U2.resize( rows - drop_rows , cols - drop_cols );
      S.resize( rows - drop_rows );
      
      int r = 0;
      for (int i=0; i<rows; i++)
	{
	  if ( row_included[i] )
	    {
	      if ( drop_cols )
		{
		  int q=0;
		  for (int j=0; j<cols; j++)
		    {
		      if ( col_included[j] )
			{
			  U2(r,q) = U(i,j);
			  ++q;
			}
		    }
		}
	      else		
		U2.row(r) = U.row(i);

	      // add stage in
	      S[r] = st2[i];
	      ++r;
	    }
	} // next row
    }

  
  const int kept = rows - blanked;
  const int nstg = stgs.size();
  const int nstg2 = stgs2.size();
  
  //
  // if leaves less than two good classes, just bail, change nothing
  //

  if ( U2.cols() < 2 ) return 0;

  if ( U2.rows() < 10 ) return 0;

  //
  // LDA-based SOAP
  //

  const int rows2 = U2.rows();

  std::cout << " S.dim = " << S.size() << " " << U2.rows() <<" " << U2.cols() << "\n";
  
  lda_t lda( S  , U2 );     

  // note: second st param means to set priors based on the full/original st[] rather
  // than the subset of unambiguous values
  
  lda_model = lda.fit( suds_t::flat_priors , st );       

  std::cout << "grops means = \n" << lda_model.means << "\n\n";
  std::cout << "grops scaling = \n" << lda_model.scaling << "\n\n";
  
  if ( ! lda_model.valid ) return 0;
  
  posteriors_t prediction = posteriors_t( lda_t::predict( lda_model , U2 ) ) ; 

  //
  // Output
  //

  if ( 0 )
    {
      int r1 = 0;
      for (int i=0; i<rows; i++)
	{
	  std::cout << (*st)[i] ;
	  
	  for (int j=0; j<cols; j++)
	    std::cout << " " << (*pp)(i,j);
	  
	  std::cout << " -> ";
	  
	  if ( row_included[i] )
	    {
	      std::cout << " " << prediction.cl[r1]
			<< " " << ( (*st)[i] != prediction.cl[r1] ? "X" : "." );
	      
	      for (int j=0; j<prediction.pp.cols(); j++)
		std::cout << " " << prediction.pp(r1,j);
	      
	      ++r1;
	    }
	  else
	    std::cout << " < -- NA -- > ";
	  
	  std::cout << "\n";
	}
    }
  
  //
  // Update (optionally, splicing back in prior epochs & posteriors)
  //
  
  int nchanged = 0;

  int r = 0;

  for (int i=0; i<rows; i++)
    {

      if ( row_included[i] )
	{

	  // replace if most-likely call has changed
	  if ( (*st)[i] != prediction.cl[r] )
	    {
	      ++nchanged;
	      (*st)[i] = prediction.cl[r] ;
	    }
	  
	  // copy posteriors back... but if cols dropped, set to 0
	  int q=0;
	  for (int j=0; j<cols; j++)
	    {
	      if ( col_included[j] )
		{
		  (*pp)(i,j) = prediction.pp(r,q);
		  ++q;
		}
	    }
	  
	  ++r;
	}
    }
  

  
  //
  // only output this in global mode
  //

  if ( 1 || global_mode )
    logger << " nstg, kep, blanked = " 
	   << nstg << " " << nstg2 << " " 
	   << kept << " " << blanked << " (tot " << ( kept + blanked ) << "\n";


  
  
  //
  // All done
  //
  
  return nchanged;
}

