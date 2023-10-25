
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

#include "models/predict.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "db/db.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "timeline/cache.h"


extern writer_t writer;
extern logger_t logger;

// init. static members
Eigen::MatrixXd model_knn_t::X = Eigen::MatrixXd::Zero( 0 , 0 );
int model_knn_t::k = 10;

prediction_t::prediction_t( edf_t & edf , param_t & param )
{

  //
  // store EDF ID for variable subst
  //

  id = edf.id;
  
  //
  // get model
  //
  
  const std::string model_file = param.requires( "model" );
  
  model.read( model_file , id );
  
  model.populate();
  
  // requires an intercept and at least one term

  if ( model.specials.find( "intercept" ) == models.specials.end() ) 
    Helper::halt( "no intercept specified in model" );
  

  // verbose output
  
  if ( param.has( "dump-model" ) )
    model.dump();
  
  //
  // always requires a cache
  //
  
  const std::string cache_name = param.requires( "cache" );
  
  if ( ! edf.timeline.cache.has_num( cache_name ) )
    Helper::halt( "cache not found: " + cache_name );
  
  cache_t<double> * cache = edf.timeline.cache.find_num( cache_name );

  //
  // kNN missing imputation
  //
  
  // allow command line to over-ride model value: 
  const std::string knn_data = param.has( "data" ) ? param.value( "data" ) : model.specials_str[ "data" ];
  
  // 0 if not defined
  int knn_n = param.has( "knn" ) ? param.requires_int( "knn" ) : model.specials[ "knn" ]; 
  
  const bool has_training_data = knn_data != "" && knn_n != 0;
  
  if ( has_training_data )
    {
      knn.load( knn_data , model.header() );
      knn.set_k( knn_n );
    }

  //
  // Z (abs) threshold to set to missing (and re-impute) 
  //
  
  const double imp_th = param.has( "th" ) ? param.requires_dbl( "th" ) : model.specials[ "th" ] ;
  const bool do_reimputation = imp_th > 0.01; 
    
  //
  // allow dropping of terms
  //

  std::set<std::string> dropped;
  if ( param.has( "drop" ) )
    dropped = param.strset( "drop" );
  
  //
  // allocate space
  //

  const int nt = model.size();

  X = Eigen::VectorXd::Zero( nt );
  Z = Eigen::VectorXd::Zero( nt );

  
  //
  // pull features from the cache
  //

  int i = 0 ;

  bool okay = true;

  // track any missing data (i.e. not present in the cache)
  
  missing.resize( nt , false );  
  int n_obs = 0;
  
  std::set<model_term_t>::const_iterator tt = model.terms.begin();
  while ( tt != model.terms.end() )
    {

      // dropped?
      if ( dropped.find( tt->label ) != dropped.end() )
	{
	  logger << "  dropping " << tt->label << "\n";
	  // this can over-ride any requirement
	  missing[i] = true;
	  ++i;
	  ++tt;
	  continue;
	}
      
      // if this a value (non-cache) term?
      if ( tt->has_value )
	{
	  double x;
	  bool valid = Helper::str2dbl( tt->value , &x ) ;
	  
	  if ( tt->value == "." || tt->value == "" || ! valid )
	    {
	      
	      logger << "  *** non-numeric/missing value specified for " << tt->label << "\n";
	      
	      if ( tt->required )
		{
		  okay = false;
		  break;
		}
	      else 
		missing[i] = true;
	    }
	  else
	    {
	      X[i] = x;
	      ++n_obs;
	    }
	  
	  // go to next term
	  ++i;
	  ++tt;
	  continue;
	}

      //
      // pull from the cache:
      //
      
      // no channels specified
      if ( tt->chs.size() == 0 && tt->pairs.size() == 0 )
	{
	  double x1;
	  if ( ! cache->fetch1( tt->cmd , tt->var , tt->strata , &x1 ) )
	    {
	      logger << "  *** could not find "
		     << tt->label << " : " << tt->cmd << " "
		     << tt->var << " "
		     << Helper::ezipam( tt->strata ) << "\n";
	      
	      // was this feature required to be non-missing?	      
	      if ( tt->required )
		{
		  // if so, a fatal error
		  okay = false;
		  break;
		}
	      else // just track and we can impute
		{
		  missing[i] = true;
		}	    
	    }

	  // add this feature 
	  X[i] = x1;
	  ++n_obs;
	  
	}
      else 
	{
	  
	  //
	  // combine features across channels (or channel pairs)
	  //
	  
	  std::vector<double> xx;

	  if ( tt->chs.size() != 0 ) // 1+ single channel specified
	    {
	      
	      // *either* look
	      // take the mean
	      
	      std::vector<std::string>::const_iterator cc = tt->chs.begin();
	      while ( cc != tt->chs.end() )
		{
		  
		  // if single channel feature, add if we have 	      
		  // single channel analysis: add CH into the strata
		  
		  std::map<std::string,std::string> ss1 = tt->strata;
		  ss1[ "CH" ] = *cc ;
		  double x1;
		  
		  // only add for this channel if present
		  if ( ! cache->fetch1( tt->cmd , tt->var , ss1 , &x1 ) )
		    {
		      logger << "  *** could not find "
			     << tt->label << " : " << tt->cmd << " "
			     << tt->var << " "
			     << Helper::ezipam( ss1 ) << "\n";
		      
		    }
		  else // rack up this channel for this feature
		    xx.push_back( x1 );
		  
		  // next channel
		  ++cc;
		}
	      
	    }
	  else 
	    {
	      
	      // CH1, CH2 scenario, specified via chs=A+B,C+D, etc
	      //   i.e. here average A+B statistic w/  C+D  statistic, etc
	      
	      // *either* look
	      // take the mean
	      
	      std::vector<std::string>::const_iterator cc = tt->pairs.begin();
	      while ( cc != tt->pairs.end() )
		{
		  
		  // add CH1 and CH2 into the strata
		  
		  std::vector<std::string> c1c2 = Helper::parse( *cc , "+" );
		  if ( c1c2.size() != 2 ) 
		    {
		      std::cout << " term " << tt->label << "\n"
				<< " term = [" << *cc << "]\n";
		      Helper::halt( "bad format for CHS=A+B,C+D,E+F" );
		    }
		  
		  std::map<std::string,std::string> ss1 = tt->strata;
		  
		  ss1[ "CH1" ] = c1c2[0];
		  ss1[ "CH2" ] = c1c2[1];
		  
		  bool retrieved = false;
		  
		  // slot
		  double x1;
		  
		  retrieved = cache->fetch1( tt->cmd , tt->var , ss1 , &x1 );
		  
		  // try swapping?
		  if ( ! retrieved ) 
		    {
		      std::map<std::string,std::string> ss2 = tt->strata;
		      ss2[ "CH2" ] = c1c2[0];
		      ss2[ "CH1" ] = c1c2[1];
		      retrieved = cache->fetch1( tt->cmd , tt->var , ss2 , &x1 );

		      // for a directed metric, need to swap sign
		      if ( tt-> directed ) x1 = - x1;
		    }
		  
		  // only add for this channel if present
		  if ( ! retrieved )
		    {
		      logger << "  *** could not find "
			     << tt->label << " : " << tt->cmd << " "
			     << tt->var << " "
			     << Helper::ezipam( ss1 ) << "\n";		      
		    }		  
		  else // rack up this channel for this feature
		    {
		      std::cout << " adding " << x1 << " for " << *cc << "\n";
		      xx.push_back( x1 );
		    }

		  // next channel pair
		  ++cc;
		}

	    }
	  

	  //
	  // check we have at least one channel (or channel pair)
	  //
	  
	  if ( xx.size() == 0 )
	    {
	      logger << "  *** could not find (for any channels) "
		     << tt->label << " : " << tt->cmd << " "
		     << tt->var << " "
		     << Helper::ezipam( tt->strata ) << "\n";

	      // fatal error?
	      if ( tt->required )
		{
		  okay = false;
		  break;
		}
	      else // just track it is missing
		{
		  missing[i] = true;
		}	      
	    }
	  else // we can add the mean across channels/pairs
	    {
	      X[i] = MiscMath::mean( xx );
	      ++n_obs;
	    }
	}

      // next term
      ++i;
      ++tt;
    }
  

  //
  // check non-missing data requirements
  //
  
  if ( n_obs < model.specials[ "minf" ] || n_obs == 0 )
    {
      logger << "  *** found " << n_obs << " non-missing features";
      if ( model.specials[ "minf" ] > 0 ) logger << " but require " << model.specials[ "minf" ] << "\n";
      okay = false;	
    }

  writer.value( "NF" , nt );
  writer.value( "NF_OBS" , n_obs );

  if ( okay )
    {
      if ( n_obs < nt )
	{
	  if ( ! knn.populated() )
	    {
	      okay = false;
	      logger << "  *** missing values, but no attached dataset for kNN imputation\n";
	    }
	}
    }


  //
  // fatality?
  //
  
  if ( ! okay )
    {
      logger << "  *** could not satisfy non-missing feature requirements... bailing\n";
      writer.value( "OKAY" , 0 );
      
      return;
    }
  writer.value( "OKAY" , 1 );
  
  
  //
  // Some checks
  //

  if ( model.mean.size() != nt )
    Helper::halt( "problem, only have " + Helper::int2str( (int)model.mean.size() ) + " means" );

  if ( model.sd.size() != nt )
    Helper::halt( "problem, only have " + Helper::int2str( (int)model.sd.size() ) + " sds" );
  
  if ( model.coef.size() != nt )
    Helper::halt( "problem, only have " + Helper::int2str( (int)model.coef.size() ) + " coefs" );
  




  
  //
  // Normalization of metrics
  //

  Z = X - model.mean ;
  
  Z = Z.array() / model.sd.array() ; 


  
  //
  // log-transformation of normalized features
  //

  bool all_logged = model.specials[ "log1p" ] ;
  
  int idx = 0,  n_tr = 0;

  tt = model.terms.begin();
  while ( tt != model.terms.end() )
    {
      if ( ! missing[idx] )
	{
	  if ( all_logged || tt->log_transform )
	    {
	      Z[idx] = ( Z[idx] > 0 ? 1 : -1 ) * log1p( fabs( Z[idx] ) );
	      ++n_tr;
	    }
	}      
      ++idx;
      ++tt;
    }

  if ( n_tr > 0 ) 
    logger << "  log1p() transformed " << n_tr << " normalized features\n";


  //
  // Missing data imputation (on Z scale)
  //
  
  if ( n_obs < nt )
    {
      logger << "  imputing missing values for " << nt - n_obs << " of " << nt << " features\n";      
      Z = knn.impute( Z , missing );
    }
  

  
  //
  // Check means, if reference data are present - assumes they are standardized
  //
  
  if ( knn.populated() )
    {      
      // original distances
      D = knn.distance( Z );
      
      // re-impute missing/weird values?
      if ( do_reimputation )
	{
	  missing2.resize( nt , false );
	  int bad = 0;
	  for (int i=0; i<nt; i++)
	    if ( fabs( D[i] ) > imp_th )
	      {
		++bad; 
		missing2[i] = true;
	      }

	  if ( bad > 0 ) 
	    {
	      logger << "  attempting re-imputation for " << bad << " features\n";
	      if ( nt - bad <  model.specials[ "minf" ] )
		{
		  logger << "  *** would imply fewer than " << model.specials[ "minf" ] << " original features remaining, bailing\n";		  
		  writer.value( "OKAY" , 0 );
		  return;
		}

	      // impute
	      Z = knn.impute( Z , missing2 );
	      
	    }
	}
            
    }

  
  
  //
  // Primary prediction
  //
  
  y = y1 = (Z.transpose() * model.coef) + model.specials[ "intercept" ]; 
       
  
  //
  // bias corrected
  //

  const bool apply_bias_correction = model.specials.find( "bias_correction_term" ) != model.specials.end() ;

  if ( apply_bias_correction ) 
    {
      if (  model.specials.find( "bias_correction_slope" ) == model.specials.end() )
	Helper::halt( "need to specify bias_correction_slope special variable\n" );
      
      if (  model.specials.find( "bias_correction_intercept" ) == model.specials.end() )
	Helper::halt( "need to specify bias_correction_intercept special variable\n" );
      
      double b = model.specials["bias_correction_slope"];
      double c = model.specials["bias_correction_intercept"];
      double x = model.specials["bias_correction_term"];
      
      y1 = y - ( b * x + c );
    }

  //
  // softplus
  //

  const bool apply_softplus = model.specials[ "softplus" ] > 0 ;
  if ( apply_softplus )
    {
      logger << "  applying softplus scaling to predicted values\n";
      
      y = log1p( exp(-fabs(y) ) ) + ( y > 0 ? y : 0 ) ; 

      if  ( apply_bias_correction ) 
	y1 = log1p( exp(-fabs(y1) ) ) + ( y1 > 0 ? y1 : 0 ) ; 
      
    }
  
  
  //
  // Primary outputs
  //

  // prediction
  logger << "\n  predicted value (Y) = " << y << "\n";
  writer.value( "Y" , y );

  // bias-corrected prediction, of model supplied
  if ( apply_bias_correction )
    {
      logger << "  bias-corrected predicted value (Y1) = " << y1 << "\n";
      writer.value( "Y1" , y1 );
    }
  
  // observed, if supplied
  if ( model.specials.find( "observed" ) != model.specials.end() )
    {
      writer.value( "YOBS" , model.specials["observed"] );
      logger << "  observed value (YOBS) = " << model.specials["observed"] << "\n";
    }
  else if ( model.specials.find( "bias_correction_term" ) != model.specials.end() )
    {
      logger << "  observed value (YOBS) = " << model.specials["bias_correction_term"] << "\n";
      writer.value( "YOBS" , model.specials["bias_correction_term"] );                                                                             
    }
  
    
  // feature level output (by FTR)
  
  output();

  //
  // all done
  //
  
}

void prediction_t::output() const
{

  int nt = model.size();

  int i = 0;  
  std::set<model_term_t>::const_iterator tt = model.terms.begin();
  while ( tt != model.terms.end() )
    {
      writer.level( tt->label , "FTR" );

      // only output if non-missing raw value
      if ( ! missing[i] ) 
	writer.value( "X" , X[i] );

      // if here, Z would have been imputed, so okay to
      // output either way
      writer.value( "Z" , Z[i] );
      
      // was a KNN run? were any features imputed
      if ( knn.populated() )
       	{
       	  // D only makes sense if non-missing
       	  if ( ! missing[i] ) 
       	    writer.value( "D" , D[i] );      
	  
       	  writer.value( "IMP" , (int)missing[i] );
	  if ( missing2.size() == nt )
	    writer.value( "REIMP" , (int)missing2[i] );
       	}

      // population/model parameters, but included for refernece
      writer.value( "M" , model.mean[i] );
      writer.value( "SD" , model.sd[i] );
      writer.value( "B" , model.coef[i] );      

      ++i;
      ++tt;
    }
  writer.unlevel( "FTR" );
  
}
