
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

// analog of POPS but for (stage-stratified) between-person classification
// DV can be binary, multiclass or quantitative trait

// similar to POPS, *except*
//   -- calculate all epoch level metrics, but summarizes by annotation (i.e. duplicate columns)
//   -- can include individual level covariates
//   -- covariates and phenotypes specified via the 'vars=<file>' mechanism
//   -- can include indiv-level stats, e.g. from HYPNO
//   -- cycle-specific analyses?  stage-specific 
//   -- all metrics are 'level-1' metrics:: i.e. no temporal smoothing etc


//   -- design so that assoc just reads in a matrix, i.e. the same as 'PSD' and '--cpt'
//      i.e. then people can include whatever metrics they like
//      and stratifications, etc...
//      but everything is done on the indiv-level
//   -- i.e. will distribute a suite of 'cmd' files along w/ a model, in order that somebody can generate
//      the necessary inputs...

#include "assoc/assoc.h"
#include "param.h"

#include "helper/logger.h"
#include "db/db.h"

lgbm_t assoc_t::lgbm;

extern logger_t logger;
extern writer_t writer;

assoc_t::assoc_t( param_t & param )
{
  
  //
  // training or test mode?
  //

  const bool training_mode = param.has( "train" ) || param.has( "training" );


  //
  // misc. options
  //
  
  lgbm.qt_mode = param.yesno( "qt" ) ;
  
  allow_missing_values = param.has( "missing" ) ? param.yesno( "missing" ) : true ; 
  
  //
  // training mode
  //

  if ( training_mode )
    {
      
      //
      // specify validation individuals
      //
      
      if ( param.has( "validation" ) )
	attach_ids( param );
  
      //
      // Get main training data
      //
      
      if ( param.has( "import" ) )
	{
	  if ( ! param.has( "save" ) )
	    Helper::halt( "'import' requires a 'save' command also" );      
	  
	  // read in from long-format text, and make training and validation IDs
	  import_training( param );      

	  // any covariates?
	  attach_covariates( param );

	  // attach a phenotype? (will be saved in the binary)
	  if ( param.has( "phe" ) ) 
	    attach_phenotypes( param );
	  
	  // save as a single binary file
	  save( param );

	  // all done now
	  return;
	}
      else
	{
	  // or we are re-loading a previously imported/save matrix (pair of matrices)
	  
	  if ( ! param.has( "load" ) )
	    Helper::halt( "training mode requires either 'import/save' or 'load' " );
	  
	  // includes features, and covariates and phenotype from import
	  load( param );
	  
	  // attach additional covariates (nb. if included in import(), then no need to specify them here again)
	  attach_covariates( param );

	  // attach a phenotype? (will be saved in the binary)                                                                                                                                  
          if ( param.has( "phe" ) )
            attach_phenotypes( param );

	}
      
      //
      // ensure we have phenotypes attached
      //

      if ( train_phe.size() == 0 )
	Helper::halt( "no phenotypes attached" );
      
      //
      // train model
      //
      
      train( param );
      
      //
      // save model
      //
      
      save_model( param );
      
	
      // all done for training mode
      return;
    }

  

  //
  // Test/prediction mode
  //

  // load model (w/ includes variable list)
  load_model( param );
  
  // load or import test data (from long text format)
  if ( param.has( "load" ) )
    load_testdata( param );
  else
    import_testdata( param );
  
  // covariates
  attach_covariates( param );
  
  // do we have real values too?  
  attach_test_phenotypes( param );
  
  // make prediction
  predict( param );

  // SHAP values
  if ( param.has( "SHAP" ) || param.has( "shap" ) )
    SHAP( param );
  
  // all done
  return;
}


void assoc_t::import_training( param_t & param )
{
  // import text/long format data, allowing for missing data (NA)
  //  assume stratifiers

  std::vector<std::string> files = param.strvector( "import" );

  // variables to import?
  std::set<std::string> vars = param.strset( "vars" );

  // stratifiers?
  std::set<std::string> strats = param.strset( "factors" );
  
  // individuals to include/exclude
  std::set<std::string> incids = param.strset( "inc-ids" );
  std::set<std::string> excids = param.strset( "exc-ids" );
  const bool check_incids = incids.size() != 0 ;
  const bool check_excids = excids.size() != 0 ;
  int skipped = 0;
  
	
  //
  // pass through all files to get the total variable list
  //
  
  var2col.clear();
  std::map<std::string,int> ind2row;
  
  for (int f=0; f<files.size(); f++)
    {
      if ( ! Helper::fileExists( Helper::expand( files[f] ) ) )
	Helper::halt( "could not open " + files[f] );

      logger << "  initial scan of " << files[f] << "\n";
      
      std::ifstream IN1( files[f].c_str() , std::ios::in );

      //
      // process header
      //
      
      std::string hline;      
      std::vector<int> sslot , vslot;      
      Helper::safe_getline( IN1 , hline );
      if ( IN1.eof() )
	Helper::halt( "problem reading from " + files[f] );
      
      std::vector<std::string> hdr = Helper::parse( hline , "\t" );

      const int ncol = hdr.size();
      
      if ( hdr.size() < 2 )
	Helper::halt( "requires at least two tab-delimited columns in " + files[f] );
      
      if ( hdr[0] != "ID" )
	Helper::halt( "expecting first column hedaer to be ID in " + files[f] );
      
      // get variables, and any factors
      for (int i=1; i<hdr.size(); i++)
	{
	  if ( strats.find( hdr[i] ) != strats.end() )
	    {
	      sslot.push_back( i );
	      model_strats.insert( hdr[i] );
	    }
	  else if ( vars.size() == 0 || vars.find( hdr[i] ) != vars.end() )
	    {
	      vslot.push_back( i );
	      model_vars.insert( hdr[i] );
	    }
	}      
      
      // process each line

      while ( 1 )
	{
	  std::string line;
	  Helper::safe_getline( IN1 , line );
	  if ( IN1.eof() ) break;
	  if ( line == "" ) continue;
	  
	  std::vector<std::string> tok = Helper::parse( line , "\t" );
	  if ( tok.size() != ncol )
	    Helper::halt( "bad line in, number of fields does not match header: "
			  + files[f] + "\n" + hline + "\n" + line );

	  // add individual as row
	  const std::string & id = tok[0];

	  // skip this person?
	  if ( check_incids && incids.find( id ) == incids.end() ) { ++skipped; continue; }
	  if ( check_excids && excids.find( id ) != excids.end() ) { ++skipped; continue; }

	  // otherwise, add to the list
	  if ( ind2row.find( id ) == ind2row.end() )
	    {
	      const int row = ind2row.size();
	      ind2row[ id ] = row;
	    }

	  // add long-format variables
	  std::string var_strats;
	  for (int ss=0; ss<sslot.size(); ss++)
	    var_strats += "_" + hdr[ sslot[ss] ] + "_" + tok[ sslot[ss] ] ;
	  
	  for (int vv=0; vv<vslot.size(); vv++)
	    {
	      
	      const std::string & var_name = hdr[ vslot[vv] ] + var_strats;

	      if ( var2col.find( var_name ) == var2col.end() )
		{
		  const int col = var2col.size();
		  var2col[ var_name ] = col;
		}
	      
	    }

	  // process next row
	}
      
      // close out this file
      IN1.close();
      
      // import next file
    }

  //
  // summarize 
  //

  const int nc = var2col.size();
  std::map<std::string,int> train2row, valid2row;
  
  std::map<std::string,int>::const_iterator ii = ind2row.begin();
  while ( ii != ind2row.end() )
    {
      if ( valid_set.find( ii->first ) != valid_set.end() )
	{
	  int n = valid2row.size();
	  valid2row[ ii->first ] = n;	  
	}
      
      else
	{
	  int n = train2row.size();
	  train2row[ ii->first ] = n;	  
	}
      
      ++ii;
    }

  const int ni_train = train2row.size();
  const int ni_valid = valid2row.size();

  train_ids.resize( ni_train );
  valid_ids.resize( ni_valid );
  std::map<std::string,int>::const_iterator kk = train2row.begin();
  while ( kk != train2row.end() )
    {
      train_ids[ kk->second ] = kk->first ;
      ++kk;
    }
  
  kk = valid2row.begin();
  while ( kk != valid2row.end() )
    {
      valid_ids[ kk->second ] = kk->first ;
      ++kk;
    }

  logger << "  expecting " << nc << " features on " << ni_train << " training and " << ni_valid << " validation observations\n";
  if ( skipped != 0 ) logger << "  skipped " << skipped << " observations due to inc-ids/exc-ids\n";

  //  req_vars.clear();
  varlist.resize( nc );
  
  std::map<std::string,int>::const_iterator vv = var2col.begin();
  while ( vv != var2col.end() )
    {
      //req_vars.insert( vv->first  );
      varlist[ vv->second ] = vv->first ;
      //std::cout << " var " << vv->first << "  --> " << vv->second << "\n"; 
      ++vv;
    }
  
  Xtrain = Eigen::MatrixXd::Constant( ni_train, nc,  NaN_value );
  Xvalid = Eigen::MatrixXd::Constant( ni_valid, nc,  NaN_value );

  
  //
  // pass through data second time, populating the matrix
  //
  
  for (int f=0; f<files.size(); f++)
    {

      logger << "  importing values from " << files[f] << "\n";
      
      std::ifstream IN1( files[f].c_str() , std::ios::in );

      // header
      
      std::string hline;      
      std::vector<int> sslot , vslot;      
      Helper::safe_getline( IN1 , hline );
      if ( IN1.eof() )
	Helper::halt( "problem reading from " + files[f] );
      std::vector<std::string> hdr = Helper::parse( hline , "\t" );
      const int ncol = hdr.size();      
      for (int i=1; i<ncol; i++)
	{
	  if ( strats.find( hdr[i] ) != strats.end() )
	    sslot.push_back( i );
	  else if ( vars.size() == 0 || vars.find( hdr[i] ) != vars.end() )
	    vslot.push_back( i );
	}      
      
      // process each line

      while ( 1 )
	{
	  std::string line;
	  Helper::safe_getline( IN1 , line );
	  if ( IN1.eof() ) break;
	  if ( line == "" ) continue;
	  
	  std::vector<std::string> tok = Helper::parse( line , "\t" );

	  // add/skip this indiv?
	  const std::string & id = tok[0];
	  if ( check_incids && incids.find( id ) == incids.end() ) continue;
	  if ( check_excids && excids.find( id ) != excids.end() ) continue;

	  // trainer (versus validiation dataset)
	  const bool is_trainer = valid_set.find( id ) == valid_set.end();
	  const int row = is_trainer ? train2row[ id ] : valid2row[ id ];
		  
	  // get stratifiers for these variables
	  std::string var_strats;
	  for (int ss=0; ss<sslot.size(); ss++)
	    var_strats += "_" + hdr[ sslot[ss] ] + "_" + tok[ sslot[ss] ] ;
	  
	  for (int vv=0; vv<vslot.size(); vv++)
	    {	      
	      const std::string & var_name = hdr[ vslot[vv] ] + var_strats;
	      const int slot = var2col[ var_name ] ;

	      double x ;
	      
	      if ( tok[ vslot[vv] ] != "NA" )
		{
		  if ( ! Helper::str2dbl( tok[ vslot[vv] ] , &x ) )
		    Helper::halt( "bad numeric value for " + id + " " + var_name + "\n" + line + "\n" );

		  // check this has not already been populated w/ a non-null value
		  const bool unpopulated = is_trainer ? std::isnan( Xtrain( row , slot ) ) :  std::isnan( Xvalid( row , slot ) );
		  
		  if ( ! unpopulated )
		    Helper::halt( "value already populated for " + id + " " + var_name + "\n" + line + "\n" );
		  
		  if ( is_trainer )
		    Xtrain( row , slot ) = x;
		  else
		    Xvalid( row , slot ) = x;
		}
	    }
	  
	  // process next row
	}
      
      // close out this file
      IN1.close();
      
      // import next file
    }

  
   // std::cout << " Xtrain\n" << Xtrain << "\n\n";
   // std::cout << " Xvalid\n" << Xvalid << "\n\n";
    
  
}

void assoc_t::attach_ids( param_t & param )
{
  // assume all individuals will be training dataset people, but allow some to be
  // put to one side (for validation)
  
  if ( param.has( "validation" ) )
    {
      const std::string filename = Helper::expand( param.value( "validation" ) );
      if ( ! Helper::fileExists( filename ) )
	Helper::halt( "could not find " + filename );
      std::ifstream IN1( filename.c_str() , std::ios::in );
      while ( 1 )
	{
	  std::string id;
	  IN1 >> id;
	  if ( IN1.eof() ) break;
	  if ( id == "" ) continue;	  
	  valid_set.insert( id );
	}
      IN1.close();
    }
  logger << "  read " << valid_set.size() << " IDs to be used as the validation set\n";
  
}


void assoc_t::attach_test_phenotypes( param_t & param )
{

  // known phenotypes optional for test dataset
  if ( ! param.has( "phe" ) ) return;
  
  phenotype_label = param.value( "phe" );
  int n_miss = 0;
  test_phe.resize( test_ids.size() );  
  for (int i=0; i<test_ids.size(); i++)
    if ( ! cmd_t::pull_ivar( test_ids[i] , phenotype_label , &(test_phe)[i] ) )
      {
	test_phe[i] = -999;
	n_miss++;
      }
  logger << "  attached " << phenotype_label << " for the test dataset";
  if ( n_miss ) logger << " (for " << test_phe.size() -n_miss << " of " << test_phe.size() << " individuals)";
  logger << "\n";
}




void assoc_t::attach_phenotypes( param_t & param )
{
  phenotype_label = param.requires( "phe" );

  // training sample
  train_phe.resize( train_ids.size() );  
  for (int i=0; i<train_ids.size(); i++)
    if ( ! cmd_t::pull_ivar( train_ids[i] , phenotype_label , &(train_phe)[i] ) )
      Helper::halt( "not phenotype " + phenotype_label + " found for " + train_ids[i] );
  logger << "  attached " << phenotype_label << " for the training dataset\n";
  
  // any validation data 
  valid_phe.resize( train_ids.size() );  
  for (int i=0; i<valid_ids.size(); i++)
    if ( ! cmd_t::pull_ivar( valid_ids[i] , phenotype_label , &(valid_phe)[i] ) )
      Helper::halt( "not phenotype " + phenotype_label + " found for " + valid_ids[i] );
  logger << "  attached " << phenotype_label << " for the validation dataset\n";
    
}

void assoc_t::save_varlist( const std::string & var_file )
{

  //
  // write out variable definitions too
  //
  

  std::ofstream O1( var_file.c_str() , std::ios::out );

  O1 << model_vars.size() << "\n";
  std::set<std::string>::const_iterator ii = model_vars.begin();
  while ( ii != model_vars.end() )
    {
      O1 << *ii << "\n";
      ++ii;
    }
  
  O1 << model_strats.size() << "\n";
  ii = model_strats.begin();
  while ( ii != model_strats.end() )
    {
      O1 << *ii << "\n";
      ++ii;
    }
  
  O1 << varlist.size() << "\n";
  std::vector<std::string>::const_iterator jj = varlist.begin();
  while ( jj != varlist.end() )
    {
      O1 << *jj << "\n";
      ++jj;
    }

  O1.close();

  logger << "  wrote variable label information to " << var_file << "\n";

}


void assoc_t::load_varlist( const std::string & var_file )
{
  
  std::ifstream I1( var_file.c_str() , std::ios::out );

  model_vars.clear();
  model_strats.clear();
  varlist.clear();

  int n;
  
  I1 >> n;
  for (int i=0;i<n;i++)
    {
      std::string x;
      I1 >> x;
      model_vars.insert(x);
    }

  I1 >> n;
  for (int i=0;i<n;i++)
    {
      std::string x;
      I1 >> x;
      model_strats.insert(x);
    }

  I1 >> n;
  for (int i=0;i<n;i++)
    {
      std::string x;
      I1 >> x;
      varlist.push_back(x);
    }
  
  I1.close();

  logger << "  read " << varlist.size() << " wide variable labels from " << var_file << "\n";

}



void assoc_t::save( param_t & param )
{

  const std::string filename = Helper::expand( param.value( "save" ) );

  // save varlist
  save_varlist( filename + ".vars" );

  // save data  
  std::ofstream OUT1( Helper::expand( filename ).c_str() , std::ios::binary | std::ios::out );

  const int nv = varlist.size();
  const int ntrain = train_ids.size();
  const int nvalid = valid_ids.size();

  logger << " writing binary data matrix, " << nv << " features, " << ntrain << " training and  " << nvalid << " validation observations\n";
  
  // number of trainers, validations
  Helper::bwrite( OUT1, nv ) ;
  Helper::bwrite( OUT1, ntrain ) ;
  Helper::bwrite( OUT1, nvalid ) ;

  // traininers
  for (int i=0; i<ntrain; i++)
    {
      // ID
      Helper::bwrite( OUT1 , train_ids[i] );
      // main phenotype (can be swapped in w/ something else)
      Helper::bwrite( OUT1 , train_phe[i] );
      // features (IVs, including covariates)
      for (int j=0; j<nv; j++)
	Helper::bwrite( OUT1 , Xtrain(i,j) );      
    }  

  // validation
  for (int i=0; i<nvalid; i++)
    {
      // ID
      Helper::bwrite( OUT1 , valid_ids[i] );
      // main phenotype (can be swapped in w/ something else)
      Helper::bwrite( OUT1 , valid_phe[i] );
      // features (IVs, including covariates)
      for (int j=0; j<nv; j++)
	Helper::bwrite( OUT1 , Xvalid(i,j) );      
    }  

  OUT1.close();
}



void assoc_t::load( param_t & param )
{

  const std::string filename = Helper::expand( param.value( "load" ) );
  
  // save varlist
  load_varlist( filename + ".vars" );

  // save data  
  std::ifstream IN1( Helper::expand( filename ).c_str() , std::ios::binary | std::ios::in );

  const int nv = Helper::bread_int( IN1 );
  const int ntrain = Helper::bread_int( IN1 );
  const int nvalid = Helper::bread_int( IN1 );

  Xtrain = Eigen::MatrixXd::Constant( ntrain, nv,  NaN_value );
  Xvalid = Eigen::MatrixXd::Constant( nvalid, nv,  NaN_value );
  
  train_ids.resize( ntrain );
  train_phe.resize( ntrain );
  valid_ids.resize( nvalid );
  valid_phe.resize( nvalid );

  // trainers
  for (int i=0; i<ntrain; i++)
    {
      train_ids[i] = Helper::bread_str( IN1 );
      train_phe[i] = Helper::bread_dbl( IN1 );
      for (int j=0; j<nv; j++)
	Xtrain(i,j) = Helper::bread_dbl( IN1 );
    }  
  
  // validation data
  for (int i=0; i<nvalid; i++)
    {
      valid_ids[i] = Helper::bread_str( IN1 );
      valid_phe[i] = Helper::bread_dbl( IN1 );
      for (int j=0; j<nv; j++)
	Xvalid(i,j) = Helper::bread_dbl( IN1 );
    }  

  IN1.close();

  logger << "  read " << nv << " variables on " << ntrain << " training and " << nvalid << " validation observations\n";
  
}


void assoc_t::load_testdata( param_t & param )
{
  // convenience function, i.e. if wanting to test the same test set under multiple
  // conditions, avoids the need to import() the text files (slower)
  
  const std::string filename = Helper::expand( param.value( "load" ) );
  
  // nb. the varlist will already have been populated from the training / load model

  std::ifstream IN1( Helper::expand( filename ).c_str() , std::ios::binary | std::ios::in );

  const int nv = Helper::bread_int( IN1 );
  const int ni = Helper::bread_int( IN1 );     // i.e. if import/save used on TEST samples only, 
  const int nvalid = Helper::bread_int( IN1 ); // and assume this will be 0 (i.e. not 'validation' specified)

  if ( nvalid != 0 )
    Helper::halt( "if loading test data, you should not have set any validation samples w/ the prior import/save" );
  
  X = Eigen::MatrixXd::Constant( ni, nv,  NaN_value );
  
  test_ids.resize( ni );
  test_phe.resize( ni );

  // test samples
  for (int i=0; i<ni; i++)
    {
      test_ids[i] = Helper::bread_str( IN1 );
      test_phe[i] = Helper::bread_dbl( IN1 );
      for (int j=0; j<nv; j++)
	X(i,j) = Helper::bread_dbl( IN1 );
    }  
  
  IN1.close();
  
  logger << "  read " << nv << " variables on " << ni << " test observations\n";
  
}

void assoc_t::train( param_t & param )
{

  // configuration required  
  lgbm.load_config( param.requires( "config" ) );

  // attach data
  lgbm.attach_training_matrix( Xtrain );
  lgbm.attach_training_qts( train_phe );

  if ( valid_ids.size() > 0 )
    {
      lgbm.attach_validation_matrix( Xvalid );
      lgbm.attach_validation_qts( valid_phe );  
    }

  
  // (TODO) attach weights?
  // ...

  // specify the number of iterations? (default = 100) 
  lgbm.n_iterations = param.has( "iter" ) ? param.requires_int( "iter" ) : 100 ; 

  // train model
  lgbm.create_booster( true );

    
}


void assoc_t::save_model( param_t & param )
{

  //
  // write out trained LGBM model
  //

  const std::string model_file = param.requires( "model" );
  
  lgbm.save_model( model_file );
    
  //
  // Also save varlist
  //

  save_varlist( model_file + ".vars" );
  
}


void assoc_t::load_model( param_t & param )
{

  const std::string model_file = param.requires( "model" );
  
  lgbm.load_model( model_file );

  logger << "  read LGBM model file from " << model_file << "\n";

  // Also load varlist

  load_varlist( model_file + ".vars" );

    
}

void assoc_t::attach_covariates( param_t & param )
{

  if ( ! param.has( "covar" ) ) return;

  const std::vector<std::string> tok = param.strvector( "covar" );
  const int ncov = tok.size();
  if ( ncov == 0 ) return;
  
  const bool training_mode = param.has( "train" ) || param.has( "training" );

  // for (int i=0; i<ncov; i++)
  //   if ( req_vars.find( tok[i] ) != req_vars.end() )
  //     Helper::halt( "covariate " + tok[i] + " is already specified as a IV - select a different label" );

  const int nv0 = training_mode ? Xtrain.cols() : X.cols();
  
  if ( training_mode )
    {
      // training sample
      Xtrain.conservativeResize( Eigen::NoChange , nv0 + ncov );
      
      for (int i=0; i<train_ids.size(); i++)
	for (int j=0; j<ncov; j++)
	  {
	    double x = NaN_value;
	    if ( cmd_t::pull_ivar( train_ids[i] , tok[j] , &x ) )
	      Xtrain( i , nv0 + j ) = x;
	  }
      
      
      // any validation data
      Xvalid.conservativeResize( Eigen::NoChange , nv0 + ncov );
      
      for (int i=0; i<valid_ids.size(); i++)
	for (int j=0; j<ncov; j++)
	  {
	    double x = NaN_value;
	    if ( cmd_t::pull_ivar( valid_ids[i] , tok[j] , &x ) )
	      Xvalid( i , nv0 + j ) = x;
	  }
      
    }
  else
    {

      // test sample(s) - the model will already include any covariates
      // and so no need to add.. we just need to match variable name

      for (int j=0; j<ncov; j++)
	{
	  const std::string covar_name = tok[j];
	  // find col/slot tok[j]
	  if ( var2col.find( covar_name ) == var2col.end() )
	    Helper::halt( "covariate " + covar_name + " not specified in the model" );

	  const int slot = var2col[ covar_name ];
	  
	  for (int i=0; i<test_ids.size(); i++)	  
	    {
	      double x = NaN_value;
	      if ( cmd_t::pull_ivar( test_ids[i] , tok[j] , &x ) )
		X( i , slot ) = x;
	    }
	}
    }
  
  // add to the varlist (for training mode only, otherwise, they will
  // already be in the varlist
  if ( training_mode )
    for (int j=0; j<ncov; j++)
      varlist.push_back( tok[j] );
  
  logger << "  attached " << ncov << " covariate(s): ";
  for (int j=0; j<ncov; j++)
    logger << " " << tok[j] ;
  logger << "\n";
}

void assoc_t::import_testdata( param_t & param )
{

  // check we have everything against req_vars... allow missing
  // unless assoc_t::allow_missing_values == F

  // track new IDs (in row order of X) in here:
  test_ids.clear();
  std::map<std::string,int> ind2row;
  

  // individuals to include/exclude
  std::set<std::string> incids = param.strset( "inc-ids" );
  std::set<std::string> excids = param.strset( "exc-ids" );
  const bool check_incids = incids.size() != 0 ;
  const bool check_excids = excids.size() != 0 ;
  int skipped = 0;
  
  // we have the set of variables read in already
  // first pass to get the # of people

  // get files
  std::vector<std::string> files = param.strvector( "import" );

  // these will be populated
  // std::set<std::string> model_vars;  long-format list of vars (e.g. PSD)
  // std::set<std::string> model_strats; long-format list of strats (e.g. CH B F)
  // std::vector<std::string> varlist: full, expanded wide-format nams P_CH_C3_F_0.5, etc
  
  var2col.clear();
  for (int i=0; i<varlist.size(); i++)
    var2col[ varlist[i] ] = i;
    
  //
  // pass through all files to get the list of individuals
  //  
  
  for (int f=0; f<files.size(); f++)
    {

      if ( ! Helper::fileExists( Helper::expand( files[f] ) ) )
	Helper::halt( "could not open " + files[f] );

      logger << "  scanning " << files[f] << "\n";
      
      std::ifstream IN1( files[f].c_str() , std::ios::in );

      // process header
      
      std::string hline;      
      Helper::safe_getline( IN1 , hline );
      if ( IN1.eof() )
	Helper::halt( "problem reading from " + files[f] );
      
      std::vector<std::string> hdr = Helper::parse( hline , "\t" );
      const int ncol = hdr.size();
      
      if ( hdr.size() < 2 )
	Helper::halt( "requires at least two tab-delimited columns in " + files[f] );
      
      if ( hdr[0] != "ID" )
	Helper::halt( "expecting first column header to be ID in " + files[f] );
      
      // process each line
      
      while ( 1 )
	{
	  std::string line;
	  Helper::safe_getline( IN1 , line );
	  if ( IN1.eof() ) break;
	  if ( line == "" ) continue;
	  
	  std::vector<std::string> tok = Helper::parse( line , "\t" );
	  if ( tok.size() != ncol )
	    Helper::halt( "bad line in, number of fields does not match header: "
			  + files[f] + "\n" + hline + "\n" + line );

	  // add individual as row?
	  const std::string & id = tok[0];

	  // skip this person?
          if ( check_incids && incids.find( id ) == incids.end() ) { ++skipped; continue; }
          if ( check_excids && excids.find( id ) != excids.end() ) { ++skipped; continue; }

	  // otherwise, add as a test subject
	  if ( ind2row.find( id ) == ind2row.end() )
	    {
	      const int row = ind2row.size();
	      ind2row[ id ] = row;
	    }
	  
	}
      
      // close out this file
      IN1.close();
      
      // import next file
    }
  

  //
  // summarize 
  //

  const int ni = ind2row.size();
  const int nv = varlist.size();

  // make storage, etc
  test_ids.resize( ni );
  std::map<std::string,int>::const_iterator ii = ind2row.begin();
  while ( ii != ind2row.end() )
    {
      test_ids[ ii->second ] = ii->first;
      ++ii;
    }
  
  logger << "  expecting " << nv << " features on " << ni << " test observations\n";
  if ( skipped != 0 ) 
    logger << "  skipped " << skipped << " observations due to inc-ids/exc-ids\n";
  
  X = Eigen::MatrixXd::Constant( ni, nv,  NaN_value );

  //
  // pass through data second time, populating the matrix
  //
  
  for (int f=0; f<files.size(); f++)
    {

      logger << "  importing values from " << files[f] << "\n";
      
      std::ifstream IN1( files[f].c_str() , std::ios::in );

      // header
      
      std::string hline;      
      std::vector<int> sslot , vslot;      
      Helper::safe_getline( IN1 , hline );
      if ( IN1.eof() )
	Helper::halt( "problem reading from " + files[f] );
      std::vector<std::string> hdr = Helper::parse( hline , "\t" );
      const int ncol = hdr.size();      
      for (int i=1; i<ncol; i++)
	{
	  if ( model_strats.find( hdr[i] ) != model_strats.end() )
	    sslot.push_back( i );
	  else if ( model_vars.find( hdr[i] ) != model_vars.end() )
	    vslot.push_back( i );
	}      
      
      // process each line

      while ( 1 )
	{
	  std::string line;
	  Helper::safe_getline( IN1 , line );
	  if ( IN1.eof() ) break;
	  if ( line == "" ) continue;
	  
	  std::vector<std::string> tok = Helper::parse( line , "\t" );

          // add/skip this indiv?
          const std::string & id = tok[0];
          if ( check_incids && incids.find( id ) == incids.end() ) continue;
          if ( check_excids && excids.find( id ) != excids.end() ) continue;

	  // get indiv row
	  const int row = ind2row[ id ];
	  
	  // get stratifiers for these variables
	  std::string var_strats;
	  for (int ss=0; ss<sslot.size(); ss++)
	    var_strats += "_" + hdr[ sslot[ss] ] + "_" + tok[ sslot[ss] ] ;
	  
	  for (int vv=0; vv<vslot.size(); vv++)
	    {	      
	      const std::string & var_name = hdr[ vslot[vv] ] + var_strats;

	      // found this variable?
	      if ( var2col.find( var_name ) != var2col.end() )
		{

		  const int slot = var2col[ var_name ] ;

		  double x ;
	      
		  if ( tok[ vslot[vv] ] != "NA" )
		    {
		      if ( ! Helper::str2dbl( tok[ vslot[vv] ] , &x ) )
			Helper::halt( "bad numeric value for " + id + " " + var_name + "\n" + line + "\n" );
		      
		      // check this has not already been populated w/ a non-null value
		      const bool unpopulated = std::isnan( X( row , slot ) ) ;
							   
		      
		      if ( ! unpopulated )
			Helper::halt( "value already populated for " + id + " " + var_name + "\n" + line + "\n" );
		  
		      X( row , slot ) = x;
		    }
		}
	    }
	  
	  // process next row
	}
      
      // close out this file
      IN1.close();
      
      // import next file
    }

  //  std::cout << " X\n" << X << "\n\n";

  
}

void assoc_t::predict( param_t & param )
{

  // only go up to iteration 'iter'?  ( 0 implies all )
  const int iter = param.has( "iter" ) ? param.requires_int( "iter" ) : 0 ;
  
  Eigen::MatrixXd Y = lgbm.predict( X , iter );
    
  const int n = Y.rows();

  //  std::cout << "Y.cols = " << Y.cols() << "\n";
  
  const int nc = Y.cols();
  
  if ( test_ids.size() != n )
    Helper::halt( "internal error in predict()" );
  
  //writer.id( "." , "." );

  for (int i=0;i<n;i++)
    {
      
      writer.id( test_ids[i] , "." );
      
      if ( test_phe.size() != 0 && test_phe[i] > -998 )
	writer.value( "OBS" , test_phe[i] );
      
      if ( nc == 1 ) // for binary & QTs
	writer.value( "PRD" , Y(i,0) );
      else // for multi-class outputs
	{
	  for (int j=0; j<nc; j++)
	    {
	      writer.level( j, "K" );
	      writer.value( "PRD" , Y(i,j) );
	    }
	  writer.unlevel( "K" );
	}

      // next obs
    }
  
  writer.id( "." , "." );

}



void assoc_t::SHAP( param_t & param )
{

  // only go up to iteration 'iter'?  ( 0 implies all )
  const int iter = param.has( "iter" ) ? param.requires_int( "iter" ) : 0 ;
  
  Eigen::MatrixXd S = lgbm.SHAP_values( X , iter );
  
  const int n = S.rows();
  const int nv = S.cols() - 1; // nb last col is expected value


  if ( test_ids.size() != n )
    Helper::halt( "internal error in predict()" );

  //  std::cout << " nv = " << nv << " " << varlist.size() << "\n";
    
  if ( nv != varlist.size() )
    Helper::halt( "internal error in predict(), varlist size" );
  
  writer.id( "." , "." );
  
  Eigen::VectorXd M = S.cwiseAbs().colwise().mean();

  //  std::cout << " S = " << S.rows() << " " << S.cols() << " " << M.size() << " " << nv << "\n";
  
  if ( M.size() != nv + 1 )
    Helper::halt( "internal error in SHAP" );

  for (int j=0; j<nv; j++)
    {
      writer.level( varlist[j] , "VAR" );
      writer.value( "SHAP" , M[j] );
    }
  writer.unlevel( "VAR" );

  // indiv level output
  
  for (int i=0;i<n;i++)
    {
      writer.id( test_ids[i] , "." );
      
      for (int j=0; j<nv; j++)
	{
	  writer.level( varlist[j] , "VAR" );
	  writer.value( "SHAP" , S(i,j) );
	}
      writer.unlevel( "VAR" );
    }
  
  writer.id( "." , "." );


}
  
#endif


