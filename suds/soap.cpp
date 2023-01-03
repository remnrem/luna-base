
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
// SOAP : Single Observation Accuracies & Probabilities
//

void suds_indiv_t::evaluate( edf_t & edf , param_t & param )
{
  
  // track ID (needed if caching for RESOAP)
  id = edf.id;

  // this impacts whether epochs w/ missing values are dropped or not  
  suds_t::soap_mode = 1;

  // ensure we do not call self_classify() from proc
  suds_t::self_classification = false;

  // verbose output
  bool epoch_level_output = param.has( "epoch" );

  // assume that we have manual staging ('true') 
  int n_unique_stages = proc( edf , param , true );
  
  //
  // Cache for RESOAP?
  //

  if ( suds_t::cache_target ) 
    {
      logger << "\n  caching " << id << " for a subsequent RESOAP\n";
      suds_t::cached = *this;
    }


  //
  // No observed stages?
  //

  if ( n_unique_stages < 2 )
    {
      logger << "  *** fewer than 2 non-missing stages for this individual, cannot complete SOAP\n";
      return;
    }


  //
  // fit LDA, and extract posteriors ( --> pp ) 
  //

  Eigen::MatrixXd pp;
  
  int dummy = self_classify( NULL , &pp );

  if ( dummy == 0 ) 
    {
      logger << "  *** not enough data/variability to fit LDA\n";
      return;
    }

  //
  // dump predictor matrix?
  //

  // to output stream
  if ( param.has( "feature-matrix" ) )
    {
      dump_predictor_matrix( edf , "" );
    }

  // as file
  if ( param.has( "dump-features" ) )
    {
      dump_predictor_matrix( edf , param.value( "dump-features" ) );
    }

  //
  // dump associations w/ stages
  //

  if ( param.has( "dump-stage-assocs" ) )
    {
      logger << "  dumping feature/SVD component stage associations to "
	     << param.value( "dump-stage-assocs" )  << "\n";
      dump_stage_associations( param.value( "dump-stage-assocs" ) );
    }

  
  //
  // dump components
  //

  if ( param.has( "dump-svd" ) )
    {
      logger << "  dumping SVD components to " << param.value( "dump-svd" )  << "\n";
      dump_svd( param.value( "dump-svd" ) );
    }
  
  //
  // output stage probabilities 
  //

  logger << "\n";
  
  const double epoch_sec = edf.timeline.epoch_length();

  const int ne_all = edf.timeline.num_epochs();

  const std::vector<std::string> & labels = suds_t::qda ? qda_model.labels : lda_model.labels ; 
  
  std::vector<std::string> final_pred = suds_t::max( pp , labels );

  summarize_kappa( final_pred , true );

  // context & stage-specific accuracies
  summarize_acc( final_pred );

  const int bad_epochs = summarize_stage_durations( pp , labels , ne_all , epoch_sec );
  
  if ( epoch_level_output )
    summarize_epochs( pp , labels , ne_all , edf );

  
  // transition reports?
  // summary of transitions
  
  if ( param.has( "trans" ) )
    {
      // specify in epochs counts
      // e.g. if 5-second epochs

      // requirements for a transition, specified in terms of old # of epochs
      const int req_left = param.has( "req-left" ) ? param.requires_int( "req-left" ) : 2 ; 
      const int req_right = param.has( "req-right" ) ? param.requires_int( "req-right" ) : 2 ;       
      
      // new epochs in seconds
      const double elen = param.requires_dbl( "trans" )  ; 

      const int show_left = param.has( "left" ) ? param.requires_int( "left" ) : 6 ; // e.g. 6 * 5 = 30 seconds
      const int show_right = param.has( "right" ) ? param.requires_int( "right" ) : 6 ; // e.g. 6 * 5 = 30 seconds
            
      summarize_transitions( pp , 
			     labels ,
			     show_left , req_left ,
			     show_right , req_right ,				      
			     elen , 
			     ne_all , edf , param );
    }
  
  
  //
  // Output annotations (of discordant epochs)
  //
  
  add_annots( pp , labels , ne_all , edf );
  
}





int suds_indiv_t::self_classify( std::vector<bool> * included , Eigen::MatrixXd * pp )
{

  if ( ! trainer )
    Helper::halt( "can only self-classify trainers (those w/ observed staging" );
  
  // assume putative 'y' and 'U' will have been constructed, and 'nve' set
  // i.e. this will be called after proc(), or from near the end of proc() 

  //
  // fit the LDA to self
  //


  fit_qlda();

  if ( suds_t::qda && ! qda_model.valid )
    return 0;
  
  if ( (!suds_t::qda) && ! lda_model.valid )
    return 0;

  
  //
  // get predictions
  //

  posteriors_t prediction;
  if ( suds_t::qda )
    prediction = posteriors_t( qda_t::predict( qda_model , U ) ) ; 
  else
    prediction = posteriors_t( lda_t::predict( lda_model , U ) ) ; 

  // save posteriors?
  if ( pp != NULL ) *pp = prediction.pp ;


  //
  // In SOAP mode, all done (we only needed the PP)
  //
  
  if ( suds_t::soap_mode || included == NULL )
    return 1;  // SOAP only cares about a non-zero return value

  //
  // Get kappa 
  //

  double kappa = MiscMath::kappa( prediction.cl , y , suds_t::str( SUDS_UNKNOWN )  );

  included->resize( nve , false );

  
  //
  // Optionally, ask whether trainer passes self-classification kappa threshold.  If not
  // make all epochs 'bad', i.e. so that this trainer will not be used
  //
  
  if ( suds_t::self_classification_kappa <= 1 )
    {
      if ( kappa < suds_t::self_classification_kappa )
	{
	  logger << "  trainer does not meet SOAP kappa " << kappa << " < " << suds_t::self_classification_kappa << "\n";
	  return 0;  // all 'included' false at this point
	}      
    }

  
  //
  // Determine 'bad' epochs
  //

  int okay = 0;

  // hard calls

  if ( suds_t::self_classification_prob > 1 )
    {
      for (int i=0;i<nve;i++)
	{
	  (*included)[i] = prediction.cl[i] == y[i] ; 
	  if ( (*included)[i] ) ++okay;
	}
    }
  else
    {
      logger << "  using threshold of PP > " << suds_t::self_classification_prob << "\n";

      // map labels to slots in PP matrix (this might be non-standard, e.g. no REM) for a
      // given trainer, and so we cannot assume canonical slot positions
      
      std::vector<std::string> labels = suds_t::qda ? qda_model.labels : lda_model.labels;
      
      std::map<std::string,int> label2slot;
      for (int j=0;j< labels.size();j++)
	label2slot[ labels[j] ] = j ;
      
      // check PP for the observated stage
      for (int i=0;i<nve;i++)
	{
	  std::map<std::string,int>::const_iterator ii = label2slot.find( y[i] );
	  if ( ii == label2slot.end() )
	    Helper::halt( "internal error in suds_indiv_t::self_classify() , unrecognized label" );
	  
	  //	  std::cout << " i " << i << "/" << nve << "  " << prediction.pp( i , ii->second ) << " " << suds_t::self_classification_prob << "\n";
	  
	  // must match hard call::: 
	  (*included)[i] = prediction.cl[i] == y[i] ;
	  
	  if ( (*included)[i] ) 
	    ++okay;
	  else // but also, must be above threshold
	    {
	      if ( prediction.pp( i , ii->second ) >= suds_t::self_classification_prob )
		{
		  (*included)[i] = true;
		  ++okay;
		}
	      else
		(*included)[i] = false;
	    }
	}
    }

  return okay;
}




void suds_indiv_t::summarize_transitions( const Eigen::MatrixXd & pp , // posterior probabilities					  					 
					  const std::vector<std::string> & labels, // column labels
					  const int show_left, const int req_left ,
					  const int show_right, const int req_right ,
					  const double elen , // new epoch length 
					  const int ne_all ,
					  edf_t & edf , param_t & param ) // total number of epochs in EDF
{
  
  // require prior staging
  const bool prior_staging = obs_stage.size() != 0 ;
  if ( ! prior_staging ) return;

  // new epoch size must be a factor of old size:
  //  i.e. integer # of new epochs in each old epoch
  const double epoch_sec = edf.timeline.epoch_length();

  double ratio = epoch_sec / elen;
  if ( fabs( ratio - round( ratio ) ) > 0.0001 )
    Helper::halt( "'trans' epoch length must be a factor of parent epoch size" );  
    
  // save this old self  
  suds_indiv_t old_self = *this;
  
  // now change epoch size to target  
  edf.timeline.set_epoch( elen , elen , 0 ) ;
  
  // and re-estimate PSD assuming no known staging ('false')
  // (this will also calculate PSC, but we will ignore this... add option to skip that in proc() in future)  
  suds_t::ignore_target_priors = true;
  
  // also clear this , as 'summarize_epochs()' will try to use it otherwise in output)
  obs_stage.clear();
  
  // re-process file
  int n_unique_stages = proc( edf , param , true ); 

  // true means has staging (i.e. not a 'target' in the SUDS sense, but 
  // but the suds_t::ignore_target_priors means this is ignored (i.e. we do 
  // not try to reference the staging (which presumably no longer matches the epoch 
  // duration)
  
  // now project & predict into self's prior PSC space;  i.e. use same model, but will just be
  // based on PSD estimated from differently-sized epochs

  posteriors_t new_staging = predict( old_self , suds_t::qda );
  
  // new posteriors : rows = epochs
  const Eigen::MatrixXd & npp = new_staging.pp;
  
  // std::cout << " npp.rows() = " << npp.rows() << "\n"
  // 	    << " es = " << epochs.size() << " " << old_self.epochs.size() << "\n";
  // Eigen::IOFormat fmt1( Eigen::StreamPrecision, Eigen::DontAlignCols );
  // std::cout << "PP\n" << npp.format( fmt1 ) << "\n";

  const std::vector<int> & e1 = old_self.epochs;
  const std::vector<int> & e2 = epochs;

  const int ne1 = ne_all;
  const int ne2 = edf.timeline.num_epochs();

  std::vector<int> ee1( ne1 , -1 );
  std::vector<int> ee2( ne2 , -1 );
  
  for (int i=0; i<e1.size(); i++) ee1[ e1[i] ] = i;
  for (int i=0; i<e2.size(); i++) ee2[ e2[i] ] = i; 

  // original staging
  const std::vector<suds_stage_t> & stages = old_self.obs_stage;
  
  int rr = ratio;

  if ( req_left < 1 || req_right < 1 )
    Helper::halt( "invalid req-left, req-right" );

  if ( show_left < 1 || show_right < 1 )
    Helper::halt( "invalid left, right" );

  // in original, flag points of transition  
  std::vector<bool> transitions( ne1 , false );
  for (int i=req_left-1; i<ne1-req_right; i++)
    {
      // set T if /next/ epoch(s) show a transition
      // N=5
      
      //     X T
      // 0 1 2 | 3 4 

      // not a transition?
      if ( stages[i] == stages[i+1] ) continue;

      // not a valid stage?
      if ( stages[i] == SUDS_UNKNOWN || stages[i] == SUDS_ARTIFACT ) continue;
      if ( stages[i+1] == SUDS_UNKNOWN || stages[i+1] == SUDS_ARTIFACT ) continue;
      
      // not a valid segment?
      if ( ee1[i] == -1 || ee1[i+1] == -1 ) continue;
      
      bool okay = true;

      // left 
      for (int j=1; j<req_left; j++)
	{
	  if ( ee1[ i-j ] == -1 )
	    {
	      okay = false;
	      continue;
	    }
	  
	  if ( stages[ i - j ] != stages[ i ] )
	    {
	      okay = false;
	      continue;
	    }
	  
	  if ( stages[ i - 1 - j ] == SUDS_UNKNOWN || stages[ i - 1 - j ] == SUDS_ARTIFACT )
	    {
	      okay = false;
              continue;
	    }

	}
      
      if ( ! okay ) continue;
      
      // right
      for (int j=1; j<req_right; j++)
	{
	  
	  if ( ee1[ i + 1 + j ] == -1 )
            {
              okay = false;
              continue;
            }
	  
	  if ( stages[ i + 1 + j ] != stages[ i + 1 ] )
	    {
	      okay = false;
	      continue;
	    }
	  
	  if ( stages[ i + 1 + j ] == SUDS_UNKNOWN || stages[ i + 1 + j ] == SUDS_ARTIFACT )
	    {
	      okay = false;
              continue;
	    }
	}
      
      
      if ( okay ) transitions[ i ] = true;      
      
    }

  
  //
  // show build transition tables
  //

  std::vector<std::string> ttype;
  std::vector<std::pair<int,int> > left, right;
  std::vector<int> left2, right2;

  int cnt = 0;
  for (int i=0; i<ne1; i++)
    {

      //std::cout << " XX " << i << "\t" << ee1[i] << "\t" << suds_t::str( stages[i] ) << "\t" << transitions[i] << "\n";	
      
      if ( transitions[i] )
	{
	  ttype.push_back( suds_t::str( stages[i] ) + "-" + suds_t::str( stages[i+1] ) );
	  left.push_back( std::make_pair( i - (req_left-1) , i ) );
	  right.push_back( std::make_pair( i + 1 , i + 1 + (req_right-1) ) );

	  // ee2 epochs
	  //  0 1 2 | 3 4 5 + 6 7 8 | 9 10 11
	  //        |       +       | 
	  //    0       1       2       3 
	  
	  // get smaller-sized epoch as the one just before the original transition
	  int key_epoch = i * rr + rr - 1 ; 
	  int key_left  = key_epoch - ( show_left - 1 );
	  int key_right = key_epoch + show_right;

	  left2.push_back( key_left );
	  right2.push_back( key_right );
			  
	  // std::cout << " ---> " << suds_t::str( stages[i] ) + "-" + suds_t::str( stages[i+1] )
	  // 	    << " -- " <<  i - (req_left-1) << " " <<  i
	  //  	    << " | " << i + 1 << " " <<  i + 1 + (req_right-1) << "\n";
	  
	  // std::cout << " ++ " << key_left << " " << key_epoch << " " << key_right << "\n";
	  
	  ++cnt;
	}
    }
  
  logger << "  found " << cnt << " valid transitions\n";
  

  //
  // now build the transition-offset means/summs
  //


  // transition type --> offset --> stage -> sum
  // transition type --> offset --> count
  
  std::map<std::string,std::map<int,std::map<std::string,double> > > tr_sums;
  std::map<std::string,std::map<int,double> > tr_counts;
  const int ns1 = labels.size();
  if ( ns1 != npp.cols() )
    Helper::halt("internal error in trans" );
    
  for (int i=0; i<cnt; i++)
    {
      //std::cout << " tr = " << i << "\n";
      
      // left
      int l = left2[i];
      int p = -show_left;
      for (int j=0; j<show_left; j++)
	{
	  //std::cout << "ee2[l] = " << ee2[l] << " p = " << p << " " << npp.rows() << " " << npp.cols() << " " << ns1 << "\n";
	  if ( ee2[l] != -1 )
	    {
	      for (int s=0; s<ns1; s++)
		tr_sums[ ttype[i] ][ p ][ labels[s] ] = npp( ee2[l] , s ); 
	      tr_counts[ ttype[i] ][ p ]++;
	    }
	  ++l;
	  ++p;
	}
      
      // right
      int r = right2[i];
      p = show_right;
      for (int j=0; j<show_right; j++)
	{
	  //std::cout << "ee2[r] = " << ee2[r] << " p = " << p << " " << npp.rows() << " " << npp.cols() << " " << ns1 << "\n";
	  if ( ee2[r] != -1 )
	    {
	      for (int s=0; s<ns1; s++)
		tr_sums[ ttype[i] ][ p ][ labels[s] ] = npp( ee2[r] , s ); 
	      tr_counts[ ttype[i] ][ p ]++;
	    }
	  --r;
	  --p;
	}
      
      // next transition
    }
  

  //
  // Now summarize
  //
  //  std::cout << " this far ... \n";
  std::map<std::string,std::map<int,std::map<std::string,double> > >::const_iterator tt = tr_sums.begin();
  while ( tt != tr_sums.end() )
    {
      writer.level( tt->first , "TTYPE" );
      
      const std::map<int,std::map<std::string,double> > & t2 = tt->second;
      std::map<int,std::map<std::string,double> >::const_iterator tt2 = t2.begin();
      while ( tt2 != t2.end() )
	{
	  writer.level( tt2->first , "OFFSET" );
	  double count = tr_counts[ tt->first ][ tt2->first ];
	  writer.value( "N" , count );
	  const std::map<std::string,double> & t3 = tt2->second;
	  std::map<std::string,double>::const_iterator tt3 = t3.begin();
	  while ( tt3 != t3.end() )
	    {
	      writer.level( tt3->first , globals::stage_strat );
	      writer.value( "PP" , tt3->second / count );
	      ++tt3;
	    }
	  writer.unlevel( globals::stage_strat );
	  ++tt2;
	}
      writer.unlevel( "OFFSET" );
      ++tt;
    }
  writer.unlevel( "TTYPE" );

  //
  // all done
  //

  
}


