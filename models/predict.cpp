
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
  
  std::set<model_term_t>::const_iterator tt = model.terms.begin();
  while ( tt != model.terms.end() )
    {

      // no channels specified
      if ( tt->chs.size() == 0 )
	{
	  double x1;
	  if ( ! cache->fetch1( tt->cmd , tt->var , tt->strata , &x1 ) )
	    {
	      okay = false;
	      break;
	    }	  
	  X[i] = x1;
	}
      else // 1+ channel specified
	{
	  // take the mean
	  std::vector<double> xx;
	  std::vector<std::string>::const_iterator cc = tt->chs.begin();
	  while ( cc != tt->chs.end() )
	    {
	      // add CH into the strata
	      std::map<std::string,std::string> ss1 = tt->strata;
	      ss1[ "CH" ] = *cc ;
	      double x1;
	      if ( ! cache->fetch1( tt->cmd , tt->var , ss1 , &x1 ) )
		{
		  logger << "  *** could not find " << tt->cmd << " "
			 << tt->var << " "
			 << Helper::ezipam( ss1 ) << "\n";
		  
		  okay = false;
		  break;
		}
	      xx.push_back( x1 );

	      // next channel
	      ++cc;
	    } 
	  
	  // take mean
	  X[i] = MiscMath::mean( xx );
	  
	}

      // next term
      ++i;
      ++tt;
    }
  

  //
  // requires complete data for now
  //
  
  if ( ! okay )
    {
      logger << "  could not find all variables... bailing\n";
      return;
    }


  //
  // Some checks
  //

  if ( model.mean.size() != nt )
    Helper::halt( "problem, only have " + Helper::int2str( model.mean.size() ) + " means" );

  if ( model.sd.size() != nt )
    Helper::halt( "problem, only have " + Helper::int2str( model.sd.size() ) + " sds" );
  
  if ( model.coef.size() != nt )
    Helper::halt( "problem, only have " + Helper::int2str( model.coef.size() ) + " coefs" );
  

  // std::cout << "X = " << X << "\n";
  // std::cout << "mean = " << model.mean << "\n";
  // std::cout << "sd = " << model.sd << "\n";
  // std::cout << "coef = " << model.coef << "\n";
    
  //
  // Normalization of metrics
  //

  Z = X - model.mean ;

  Z = X.array() / model.sd.array() ; 

  //
  // Output
  //

  output();

  //
  // Prediction
  //

  // raw
  y = y1 = (Z.transpose() * model.coef) + model.specials[ "model_intercept" ]; 
  logger << "  predicted value = " << y << "\n";
  writer.value( "Y" , y );
  
  // bias corrected

  if ( model.specials.find( "bias_correction_term" ) != model.specials.end() )
    {
      double b = model.specials["bias_correction_slope"];
      double c = model.specials["bias_correction_intercept"];
      double x = model.specials["bias_correction_term"];
      
      y1 = y - ( b * x + c );
      logger << "  bias-corrected predicted value = " << y1 << "\n";
      writer.value( "Y1" , y );
      writer.value( "X" , x );
    }

  
  
}


void prediction_t::output() const
{

  int nt = model.size();

  int i = 0;  
  std::set<model_term_t>::const_iterator tt = model.terms.begin();
  while ( tt != model.terms.end() )
    {
      writer.level( tt->label , "TERM" );
      writer.value( "X" , X[i] );
      writer.value( "Z" , Z[i] );
      writer.value( "M" , model.mean[i] );
      writer.value( "SD" , model.sd[i] );
      writer.value( "B" , model.coef[i] );      
      ++i;
      ++tt;
    }
  writer.unlevel( "TERM" );
  
}

