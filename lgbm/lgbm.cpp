
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


// compilation notes on ERIS
// module load cmake/3.16.0
// export CXX=/apps/lib-osver/gcc/6.3.0/bin/g++ CC=/apps/lib-osver/gcc/6.3.0/bin/gcc
//  prior to cmake/make process as outlined

#ifdef HAS_LGBM

#include "lgbm/lgbm.h"
#include "param.h"
#include "stats/eigen_ops.h"
#include <fstream>

#include "helper/logger.h"
#include "db/db.h"

extern logger_t logger;
extern writer_t writer;

void lgbm_cli_wrapper( param_t & param )
{

  // training
  //   attach configuration file
  //          training data
  //          validation data
  //          label/weight file
  //   save model

  // query model: e.g. SHAP values
  
  // predict
  //   attach model
  //          test data
  //   save posteriors/labels
  
  const bool has_training = param.has( "train" );
  const bool has_training_weights = param.has( "train-weights" );

  const bool has_validation = param.has( "valid" );
  const bool has_validation_weights = param.has( "valid-weights" );

  const bool has_label_weights = param.has( "weights" );

  if ( has_label_weights && ( has_training_weights || has_validation_weights ) )
    Helper::halt( "can only specify weights or train-weights/valid-weights" );
  
  const bool has_test = param.has( "test" );
  const bool has_config = param.has( "config" );
  const std::string model_file = param.requires( "model" );

  const bool out_shap = param.has( "SHAP" ) || param.has( "shap" );
  const bool qt_mode = param.has( "qt" );
  
  if ( has_training && has_test ) Helper::halt( "can only specify train or test" );
  if ( ! ( has_training || has_test ) ) Helper::halt( "no train or test data attached" );
  if ( has_validation && ! has_training ) Helper::halt( "can only specify valid with train" );
  
  // generic input formats:
  //  - for test, LABEL typically unknown ('.') 
  // ID LABEL F1 F2 ...

  
  lgbm_t lgbm;

  // 
  // classification or regression?
  //

  lgbm.qt_mode = qt_mode;

  //
  // attach configuration file
  //
  
  if ( has_config )
    lgbm.load_config( param.value( "config" ) );
  
  //
  // load training
  //

  if ( has_training )
    {
      lgbm.load_training_data( param.value( "train" ) );
    
      logger << "  attached training data ("
       	     << lgbm_t::rows( lgbm.training ) << " x "
       	     << lgbm_t::cols( lgbm.training ) << " ) from "
       	     << param.value( "train" ) << "\n";
    }
  
  //
  // validation data
  //

  if ( has_validation )
    {
      lgbm.load_validation_data( param.value( "valid" ) );

      logger << "  attached validation data ("
	     << lgbm_t::rows( lgbm.validation ) << " x "
	     << lgbm_t::cols( lgbm.validation ) << " ) from "
	     << param.value( "valid" ) << "\n";
    }
  
  //
  // Weights?
  //
  
  // per-label weights (applied to both datasets)

  if ( has_label_weights )
    {
      if ( qt_mode )
	Helper::halt( "cannot apply label weights in QT mode" );
      
      lgbm_label_t labels( param.value( "weights" ) );
      logger << "  applying label-weights from " << param.value( "weights" ) << "\n";
      
      if ( has_training ) 
	lgbm.add_label_weights( lgbm.training , &lgbm.training_weights, labels );
      
      if ( has_validation ) 
	lgbm.add_label_weights( lgbm.validation , &lgbm.validation_weights, labels );
    }
  
  // per-observation weight file: training
  if ( has_training_weights )
    {
      logger << "  attached training weights from " << param.has( "train-weights" ) << "\n";
      lgbm.load_weights( lgbm.training , &lgbm.training_weights, param.value( "train-weights" ) );      
    }

  // per-observation weight file: validation
  if ( has_validation_weights )
    {
      logger << "  attached validation weights from " << param.has( "valid-weights" ) << "\n";
      lgbm.load_weights( lgbm.validation , &lgbm.validation_weights, param.value( "valid-weights" ) );      
    }


  //
  // Apply weights
  //
  
  if ( has_label_weights || has_validation_weights ) 
    {
      if ( has_training ) 
	lgbm.apply_weights( lgbm.training , &lgbm.training_weights );

      if ( has_validation ) 
	lgbm.apply_weights( lgbm.validation , &lgbm.validation_weights );
      
    }

  //
  // Train and save model
  //

  if ( has_training ) 
    {
      lgbm.create_booster();
      
      lgbm.save_model( model_file );
      
      // all done
      return;
    }

  
  //
  // Prediction mode
  //
  
  const bool has_header = param.has( "header" ) ? param.yesno( "header" ) : true ;
  const bool has_ids = param.has( "ids" ) ? param.yesno( "ids" ) : true ;
  const bool has_labels = param.has( "labels" ) ? param.yesno( "labels" ) : true ;

  std::vector<std::string> headers, ids, labels;

  Eigen::MatrixXd X = eigen_ops::load_mat( param.requires( "test" ) ,
					   has_header ? &headers : NULL ,
					   has_ids ? &ids : NULL ,
					   has_labels ? &labels : NULL );
  
  logger << "  read test data (" << X.rows() << " x " << X.cols()
	 << ") from " << param.requires( "test" )  << "\n";

  //
  // Load model and predict
  //

  lgbm.load_model( model_file );
  
  Eigen::MatrixXd P = lgbm.predict( X );
  
  //
  // Output predictions
  //

  const int nobs = P.rows();
  const int ncol = P.cols();
  //  std::cout << "P\n" << P << "\n";
  
  
 
}

// https://github.com/microsoft/LightGBM/blob/master/tests/c_api_test/test_.py
// https://github.com/microsoft/LightGBM/issues/2625

bool lgbm_t::create_booster( const bool verbose )
{
  // LGBM_BoosterCreate
  // LGBM_BoosterAddValidData(booster, test)
  // loop (iterations...)
  //      LGBM_BoosterUpdateOneIter()
  //      LGBM_BoosterGetEval()
  //      LGBM_BoosterSaveModel(
  //      LGBM_BoosterFree(booster)

  int flag = LGBM_BoosterCreate( training ,
				 params.c_str() ,
				 &booster );
  
  if ( flag )
    Helper::halt( "problem creating this file" );

  has_booster = true;
  
  //
  // add validation data 
  //

  if ( has_validation )
    {  
      flag = LGBM_BoosterAddValidData( booster , validation );
      
      if ( flag )
	Helper::halt( "problem adding validation data" );
    }
  
  
  //
  // iterations
  //

  int num_validation_eval_metrics = 0;

  LGBM_BoosterGetEvalCounts( booster , &num_validation_eval_metrics );
  
  for (int i = 0; i < n_iterations ; i++)
    {
      int is_finished;
      
      int flag = LGBM_BoosterUpdateOneIter( booster , &is_finished );
      
      if ( flag )
	Helper::halt( "problem iterating training model" );
      
      if ( is_finished == 1 )
	{
	  logger << "  finished in " << i+1 << " iterations\n";
	  break;
	}
      
      //
      // Evaluation
      //
      
      int out_len;
      std::vector<double> eval( num_validation_eval_metrics , 0 );
      
      flag = LGBM_BoosterGetEval( booster ,
				  0 ,  //  dataset index == 0 training, 1 = 1st validation dataset, etc
				  &out_len ,
				  eval.data() );

      if ( flag ) 
	Helper::halt( "problem evaluating training data" );

      int out_len_valid;
      std::vector<double> eval_valid( num_validation_eval_metrics , 0 );
      
      if ( has_validation )
	{
	  flag = LGBM_BoosterGetEval( booster ,
				      1 ,  //  dataset index == 0 training, 1 = 1st validation dataset, etc
				      &out_len_valid ,
				      eval_valid.data() );

	  if ( flag ) 
	    Helper::halt( "problem evaluating validation data" );
	}
    
      logger << " iteration " << i+1 << ": training =";
      
      for (int j=0; j<out_len; j++) 
	logger << " " << eval[j];
      if ( has_validation ) 
	{
	  logger << " validation =";
	  for (int j=0; j<out_len_valid; j++)
	    logger << " " << eval_valid[j];
	}
      logger << "\n";
      
      
      // track in DB?
      if ( verbose )
	{
	  writer.level( i+1 , "ITER" );
	  for (int j=0; j<out_len; j++)
	    {
	      writer.level( j+1 , "METRIC" );
	      writer.value( "TRAINING" ,  eval[j] );
	      if ( has_validation )
		writer.value( "VALIDATION" ,  eval_valid[j] );
	    }
	  writer.unlevel( "METRIC" );
	}
      
    }

  if ( verbose )
    writer.unlevel( "ITER" );
  
  //  - LGBM_BoosterGetEvalNames first to get the names of evaluation metrics
  //  - pre-allocate memory for out_results (length by LGBM_BoosterGetEvalCounts)

  
  return true;
}

bool lgbm_t::attach_training_matrix( const Eigen::MatrixXd & X )
{
  
  int res = LGBM_DatasetCreateFromMat( X.data() , 
				       C_API_DTYPE_FLOAT64 ,
				       X.rows() ,
				       X.cols() ,
				       0 , // col-major
				       params.c_str() ,
				       nullptr ,
				       &training );

  if ( res ) Helper::halt( "problem attaching training data" );

  // set all weights to 1.0
  reset_weights( training , &training_weights );

  has_training = true;
  return true;

}

bool lgbm_t::attach_validation_matrix( const Eigen::MatrixXd & X )
{
  
  int res = LGBM_DatasetCreateFromMat( X.data() , 
				       C_API_DTYPE_FLOAT64 ,
				       X.rows() ,
				       X.cols() ,
				       0 , // col-major
				       params.c_str() ,
				       training ,
				       &validation );

  if ( res ) Helper::halt( "problem attaching validation data" );

  // set all weights to 1.0
  reset_weights( validation , &validation_weights );

  has_validation = true;
  return true;

}


bool lgbm_t::load_training_data( const std::string & f )
{
  // load training data
  std::string filename = Helper::expand( f );  

  if ( ! Helper::fileExists( filename ) )
    Helper::halt( "could not open " + filename );

  //  std::cout << "params [" << params << "]\n";
  
  int res = LGBM_DatasetCreateFromFile( filename.c_str() , 
 					params.c_str() , 
					nullptr , 
					&training );

  if ( res ) Helper::halt( "problem loading training data" );

  // set all weights to 1.0
  reset_weights( training , &training_weights );

  has_training = true;  
  return true;
}

bool lgbm_t::load_validation_data( const std::string & f )
{
  // load validation data
  //  - use training dataset as reference  
  std::string filename = Helper::expand( f );
  
  if ( ! Helper::fileExists( filename ) )
    Helper::halt( "could not open " + filename );  
  
  int res = LGBM_DatasetCreateFromFile( filename.c_str() , 
 					params.c_str() , 
					training , 
					&validation );  
  if ( res ) Helper::halt( "problem loading validation data" );

  // set all weights to 1.0
  reset_weights( validation , &validation_weights );

  has_validation = true; 
  return true;
}


bool lgbm_t::attach_training_labels( const std::vector<int> & labels )
{
  const int n = labels.size();
  std::vector<float> fl( n );
  for (int i=0; i<n; i++) fl[i] = labels[i];

  int res = LGBM_DatasetSetField( training , 
				  "label" ,
				  fl.data() ,
				  n ,
				  C_API_DTYPE_FLOAT32 );

  if ( res )
    Helper::halt( "problem attaching training labels" );

  return true;
}


bool lgbm_t::attach_training_qts( const std::vector<double> & qts )
{
  const int n = qts.size();
  std::vector<float> fl( n );
  for (int i=0; i<n; i++) fl[i] = qts[i];
  
  int res = LGBM_DatasetSetField( training , 
				  "label" ,
				  fl.data() ,
				  n ,
				  C_API_DTYPE_FLOAT32 );

  if ( res )
    Helper::halt( "problem attaching training labels" );

  return true;
}


bool lgbm_t::attach_validation_labels( const std::vector<int> & labels )
{
  const int n = labels.size();
  std::vector<float> fl( n );
  for (int i=0; i<n; i++) fl[i] = labels[i];

  int res = LGBM_DatasetSetField( validation ,
				  "label" ,
				  fl.data() ,
				  n ,
				  C_API_DTYPE_FLOAT32 );

  if ( res )
    Helper::halt( "problem attaching validation labels" );

  return true;
}

bool lgbm_t::attach_validation_qts( const std::vector<double> & qts )
{
  const int n = qts.size();
  std::vector<float> fl( n );
  for (int i=0; i<n; i++) fl[i] = qts[i];

  int res = LGBM_DatasetSetField( validation ,
				  "label" ,
				  fl.data() ,
				  n ,
				  C_API_DTYPE_FLOAT32 );

  if ( res )
    Helper::halt( "problem attaching validation labels" );

  return true;
}

bool lgbm_t::reset_weights( DatasetHandle d , std::vector<float> * w )
{
  const int n = rows(d);
  w->resize(n);
  for (int i=0; i<n; i++) (*w)[i] = 1.0;
  return true;
}

// set a weight column to a dataset (from file)
bool lgbm_t::load_weights( DatasetHandle d , std::vector<float> * w , const std::string & f )
{
  std::string filename = Helper::expand( f );
  
  if ( ! Helper::fileExists( filename ) )
    Helper::halt( "could not attach weight file " + filename );
  
  w->clear();
  std::ifstream IN1( filename.c_str() , std::ios::in );
  while ( 1 )
    {
      float x;
      IN1 >> x;
      if ( IN1.bad() || IN1.eof() ) break;
      w->push_back(x);      
    }
  IN1.close();
  
  logger << "  reading " << w->size() << " weights from " << filename << "\n";
  
  return true;
}


// load model from a file
bool lgbm_t::load_model( const std::string & f )
{
  std::string filename = Helper::expand( f );  
  if ( ! Helper::fileExists( filename ) ) Helper::halt( "could not open " + filename );
  int out_num_iterations;  
  int temp = LGBM_BoosterCreateFromModelfile( filename.c_str() ,
					      &out_num_iterations ,
					      &booster );
  has_booster = true;  
  logger << "  read model from " << filename << " (" << out_num_iterations << " iterations)\n";
  return true;
}

// load model from a string
bool lgbm_t::load_model_string( const std::string & str )
{  
  int out_num_iterations;  
  int res = LGBM_BoosterLoadModelFromString( str.c_str() ,
					     &out_num_iterations ,
					     &booster );
  if ( res ) Helper::halt( "problem in lgmb_t::load_model()" );
  logger << "  attached model (" << out_num_iterations << " iterations)\n";  
  return true;
}


bool lgbm_t::save_model( const std::string & filename )
{
  int res = LGBM_BoosterSaveModel( booster , 
 				   0 , // start_iteration – Start index of the iteration that should be saved 
 				   0 , // num_iteration – Index of the iteration that should be saved, <= 0 means save all
 				   C_API_FEATURE_IMPORTANCE_SPLIT , // or C_API_FEATURE_IMPORTANCE_GAIN ??
 				   Helper::expand( filename ).c_str() );
  
  if ( res ) Helper::halt( "problem in lgmb_t::save_model()" );  
  logger << "  saved model file to " << filename << "\n";  
  return true;
}




Eigen::MatrixXd lgbm_t::predict( const Eigen::MatrixXd & X , const int final_iter )
{
  
  // std::cout << "X dim = " << X.rows() << " " << X.cols() << "\n";
  
  // std::cout << X << "\n\n";


  if ( ! has_booster )
    Helper::halt( "no model defined" );

  const void * p = static_cast<const void*>(X.data());
  
  // results
  // std::cout << " qt_mode = " << qt_mode << "\n";
  // std::cout << " bbb " << lgbm_t::classes( booster ) << "\n";

  int num_classes = qt_mode ? 1 : lgbm_t::classes( booster );
  int num_obs = X.rows();

  //  std::cout << " num , obs = " << num_classes <<" " << num_obs<< "\n";
  
  int64_t out_len = num_classes * num_obs;

  //  std::cout << " out result(1) " << out_len << "\n";

  // note: returns row-major storage, so read into transposed matrix
  //  and transpose on return (below)
  Eigen::MatrixXd R( num_classes , num_obs );
  
  double * out_result = R.data();

  int flag = LGBM_BoosterPredictForMat( booster ,
					p ,
					C_API_DTYPE_FLOAT64 , // or C_API_DTYPE_FLOAT32 
					X.rows() ,
					X.cols() ,
					0 , // 1=row_major, 0=col-major
					C_API_PREDICT_NORMAL , //
					0 , // start_iteration
					final_iter , // Number of iteration for prediction, <= 0 means no limit
					params.c_str() ,
					&out_len ,
					out_result );
  
  // std::cout << " done\n";
  // std::cout << " out result(2) " << out_len << "\n";

  if ( flag )
    Helper::halt( "issue w/ prediction" );

  // std::cout << " num_classes = " << num_classes << "\n"
  // 	    << " qt_mode = " << qt_mode << "\n";
  
  // for binary classificaiton, make a two-col matrix
  // (i.e. same as for multiclass)
  if ( num_classes == 1 && ! qt_mode )
    {
      R.conservativeResize( 2 , Eigen::NoChange );
      for (int i=0; i<R.cols(); i++)
	R(1,i) = 1 - R(0,i);
    }

  // std::cout << "R = " << R.rows() <<" " << R.cols() <<"\n";
  
  // std::cout << " R(t) \n" << R.transpose() << "\n";
  
  return R.transpose();
  
}


Eigen::MatrixXd lgbm_t::SHAP_values( const Eigen::MatrixXd & X , const int final_iter )
{

  const void * p = static_cast<const void*>(X.data());
  
  int64_t out_len;

  int flag = LGBM_BoosterCalcNumPredict( booster , 
					 1 ,  // number of rows- i.e. just a multiplicative factor
					 C_API_PREDICT_CONTRIB , // SHAP values
					 0 ,  // start iteration
					 final_iter ,  // end (0 -> no limit)
					 &out_len );
  
  if ( flag )
    Helper::halt( "issue w/ getting SHAP values" );
  
  // for feature contributions, its length is equal to num_class * num_data * (num_feature + 1).
  int num_classes = qt_mode ? 1 : lgbm_t::classes( booster );
  int num_obs = X.rows();
  int num_features = X.cols();
  int64_t out_len2 = num_classes * num_obs * ( num_features - 1 );
  int64_t out_len3;
  
  std::vector<double> R( out_len * num_obs , 0 );
  
  double * out_result = R.data();

  flag = LGBM_BoosterPredictForMat( booster ,
				    p ,
				    C_API_DTYPE_FLOAT64 , // or C_API_DTYPE_FLOAT32
				    X.rows() ,
				    X.cols() ,
				    0 , // 1=row_major, 0=col-major
				    C_API_PREDICT_CONTRIB , // SHAP values
				    0 , // start_iteration
				    final_iter , // Number of iteration for prediction, <= 0 means no limit
				    params.c_str() ,
				    &out_len3 ,
				    out_result );

  
  if ( flag )
    Helper::halt( "issue w/ getting SHAP values" );
  
  
  // Epoch
  //  then class
  //    then features... (last col = expected value)
  
  //  std::cout << "nums " << out_len3 <<  " " << num_obs << " " << num_classes << " " <<  num_features  << "\n";
  if ( out_len3 != num_obs * num_classes * (num_features+1) )
    Helper::halt( "internal error in SHAP()" );
  
  Eigen::MatrixXd SV = Eigen::MatrixXd( num_obs , num_classes * (num_features+1) );
  
  int cnt = 0;
  for (int i=0; i<num_obs; i++)
    for (int j=0; j<num_classes; j++)
      for (int k=0; k<num_features+1; k++)
	SV(i,j*num_features+k) = R[cnt++];

  return SV;

}




std::string lgbm_t::parse_config( const std::string & f )
{
  
  std::string s;

  std::string filename = Helper::expand( f );
  
  if ( ! Helper::fileExists( filename ) )
    Helper::halt( "could not open " + filename );
  
  std::ifstream IN1( filename.c_str() , std::ios::in );
  
  while ( 1 )
    {
      std::string line;
      Helper::safe_getline( IN1 , line );
      if ( IN1.bad() || IN1.eof() ) break;
      if ( line == "" ) continue;
      if ( line[0] == '#' ) continue;

      // remove spaces
      line.erase(remove(line.begin(), line.end(), ' '), line.end());      
      s += line + " ";
    }
  IN1.close();

  return s;

}



int lgbm_t::classes( BoosterHandle b )
{
  int out = 0;
  int res = LGBM_BoosterGetNumClasses( b , &out );
  if ( res ) Helper::halt( "internal error in lgbm_t::classes()" );
  //  std::cout << " booster num_classes = " << out << "\n";
  return out;
}

int lgbm_t::cols( DatasetHandle d )
{
  int out = 0;
  int res = LGBM_DatasetGetNumFeature( d , &out );
  if ( res ) Helper::halt( "internal error in lgbm_t::cols()" );
  return out;
}

int lgbm_t::rows( DatasetHandle d )
{
  int out = 0;
  int res = LGBM_DatasetGetNumData( d , &out );
  if ( res ) Helper::halt( "internal error in lgbm_t::rows()" );
  return out;
}

int lgbm_t::label_column( DatasetHandle d )
{
  return -1;
}


std::vector<int> lgbm_t::labels( DatasetHandle d )
{
  const int n = rows(d);
  int out_len = 0;
  const void * out_ptr;
  int out_type;
  
  int res = LGBM_DatasetGetField( d ,
				  "label" , 
				  &out_len ,
				  &out_ptr ,
				  &out_type );
  
  if ( res ) Helper::halt( "problem in lgbm_t::labels" );
  if ( out_len != n ) Helper::halt( "internal error in lgbm_t::labels()" );
  
  // return labels as ints
  std::vector<int> rv(n);
  
  if ( out_type == C_API_DTYPE_FLOAT32 )
    {
      float * pp = (float*)out_ptr;      
      for (int i=0; i<n; i++) rv[i] = (int)(*(pp++));
    }
  
  if ( out_type == C_API_DTYPE_FLOAT64 )
    {
      double * pp = (double*)out_ptr;
      for (int i=0; i<n; i++) rv[i] = (int)(*(pp++));
    }

  if ( out_type == C_API_DTYPE_INT32 )
    {
      int * pp = (int*)out_ptr;
      for (int i=0; i<n; i++) rv[i] = *(pp++);      
    }

  return rv;
}


std::vector<double> lgbm_t::qts( DatasetHandle d )
{
  const int n = rows(d);
  int out_len = 0;
  const void * out_ptr;
  int out_type;
  
  int res = LGBM_DatasetGetField( d ,
				  "label" , 
				  &out_len ,
				  &out_ptr ,
				  &out_type );
  
  if ( res ) Helper::halt( "problem in lgbm_t::labels" );
  if ( out_len != n ) Helper::halt( "internal error in lgbm_t::labels()" );
  
  // return labels as doubles
  std::vector<double> rv(n);
  
  if ( out_type == C_API_DTYPE_FLOAT32 )
    {
      float * pp = (float*)out_ptr;      
      for (int i=0; i<n; i++) rv[i] = (double)(*(pp++));
    }
  
  if ( out_type == C_API_DTYPE_FLOAT64 )
    {
      double * pp = (double*)out_ptr;
      for (int i=0; i<n; i++) rv[i] = (double)(*(pp++));
    }

  if ( out_type == C_API_DTYPE_INT32 )
    {
      int * pp = (int*)out_ptr;
      for (int i=0; i<n; i++) rv[i] = (double)(*(pp++));
    }

  return rv;
}


std::vector<double> lgbm_t::weights( DatasetHandle d )
{
  const int n = rows(d);
  int out_len = 0;
  const void * out_ptr;
  int out_type;
  
  int res = LGBM_DatasetGetField( d ,
				  "weight" , 
				  &out_len ,
				  &out_ptr ,
				  &out_type );
  
  if ( res ) Helper::halt( "problem in lgbm_t::labels" );
  if ( out_len != n ) Helper::halt( "internal error in lgbm_t::labels()" );
  
  // return labels as doubles
  std::vector<double> rv(n);
  
  if ( out_type == C_API_DTYPE_FLOAT32 )
    {
      float * pp = (float*)out_ptr;      
      for (int i=0; i<n; i++) rv[i] = *(pp++);
    }
  
  if ( out_type == C_API_DTYPE_FLOAT64 )
    {
      double * pp = (double*)out_ptr;
      for (int i=0; i<n; i++) rv[i] = *(pp++);
    }

  if ( out_type == C_API_DTYPE_INT32 )
    {
      int * pp = (int*)out_ptr;
      for (int i=0; i<n; i++) rv[i] = *(pp++);      
    }

  return rv;
}


std::vector<std::string> lgbm_t::features( DatasetHandle d )
{
  std::vector<std::string> f;
  return f;
}


bool lgbm_t::apply_weights( DatasetHandle d , std::vector<float> * w )
{
  
  // apply weights
  int res = LGBM_DatasetSetField( d ,
				  "weight" , 
				  w->data() ,
				  w->size() , 
				  C_API_DTYPE_FLOAT32 );
  
  if ( res )
    Helper::halt( "problem attaching weights" );
  
  return true;
  
}

bool lgbm_t::add_label_weights( DatasetHandle d , std::vector<float> * w , const lgbm_label_t & l )
{
  
  // get labels;
  std::vector<int> lab = lgbm_t::labels( d );
  const int n = lgbm_t::rows( d );
  
  for (int i=0; i<n; i++)
    {
      if ( lab[i] < 0 || lab[i] >= l.n )
	Helper::halt( "internal error in lgbm_t::apply_label_weights() " );
      
      // multiplicative weight update
      (*w)[i] *= l.weight[ lab[i] ];
    }
  
  return true;
}

bool lgbm_t::add_block_weights( DatasetHandle d , 
				std::vector<float> * w, 
				const std::vector<uint64_t> & Istart , 
				const std::map<uint64_t,float > & wtable )
{

  const int nind = Istart.size();  
  if ( nind == 0 ) return false;
  
  const int ndat = rows(d);
  const int nrow = w->size();
  if ( ndat != nrow ) 
    Helper::halt( "internal problem in lgbm_t::add_block_weights()");
  
  // go up to the N-1 person
  for (int i=0; i<nind-1; i++)
    {
      //std::cout <<" looking for " << Istart[ i ] << "\n";
      std::map<uint64_t,float>::const_iterator ff = wtable.find( Istart[ i ] ) ;
      if ( ff != wtable.end() )
	{
	  int s1 = Istart[i];
	  int s2 = Istart[i+1];
	  //std::cout << " setting " << Istart[i] <<" " << s1 << " to " << s2 <<" " << ff->second << "\n";
	  for (int j=s1; j<s2; j++)
	    (*w)[j] *= ff->second;
	}
    }
  
  // go from last person to the end (nb. if case of only one person, this will
  // also be the first person
  std::map<uint64_t,float>::const_iterator ff = wtable.find( Istart[ nind -1  ] ) ;
  if ( ff != wtable.end() )
    for (int j= Istart[ nind-1 ] ; j<nrow; j++)
      (*w)[j] *= ff->second;
  
  return true;
}

void lgbm_t::load_pops_default_config()
{

  params = 
    " boosting_type = gbdt"
    " objective = multiclass"
    " metric = multi_logloss"
    " num_class = 5"
    " metric_freq = 1"
    " is_training_metric = true"
    " max_bin = 255"
    " early_stopping = 10"
    " num_trees = 100"
    " learning_rate = 0.05"
    " num_leaves = 31";  
  
}

   
#endif





