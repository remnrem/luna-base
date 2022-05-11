
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

#include "assoc/massoc.h"
#include "db/db.h"
#include "helper/logger.h"
#include "stats/eigen_ops.h"
#include "stats/statistics.h"
#include <fstream>

lgbm_t massoc_t::lgbm;

extern logger_t logger;
extern writer_t writer;


//
// when called externally, i.e. for testing/prediction
//

massoc_t::massoc_t( param_t & param )
{

  // runs in 1 of 6 modes
  //   train
  //   split
  //   merge
  //   dump
  //   rows
  //   test (default) [ also, this is 'dump' mode ] 

  const bool split_mode = param.has( "split" );
  const bool merge_mode = param.has( "merge" );
  const bool train_mode = param.has( "train" );
  const bool rows_mode = param.has( "rows" );
  const bool dump_mode = param.has( "dump" );

  // the default:
  const bool test_mode = param.has( "test" ) || ! ( split_mode || train_mode || merge_mode || dump_mode || rows_mode );
  
  if ( train_mode ) mode = 1;
  else if ( test_mode ) mode = 2;
  else mode = 0;

  if ( test_mode + train_mode + split_mode + merge_mode + dump_mode + rows_mode > 1 )
    Helper::halt( "can only specify one of split, merge, dump, train or test" );
  
  if ( train_mode && ! param.has( "phe" ) )
    Helper::halt( "'phe' required for 'train' mode" );
  
  //
  // dump mode - take a single file and output row IDs
  //

  if ( rows_mode )
    {
      const std::string infile = Helper::expand( param.requires( "load" ) );

      load( infile , 1 );

      const int nrow = training_iids.size();

      std::cout << "ID\n";
      for (int i=0; i<nrow; i++)
	std::cout << training_iids[i] << "_" << training_ids[i] << "_" << training_eids[i] << "\n";

      return;
      
    }
  
  //
  // merge mode - take two files and combine 
  //

  if ( merge_mode )
    {
      std::string infile1 = param.requires( "load1" );
      std::string infile2 = param.requires( "load2" );
      std::string outfile = param.requires( "save" );
      
      logger << "  merging " << infile1 << " and " << infile2 << " -> " << outfile << "\n";
      
      //  1 --> means split mode, so load all into Xtrain
      load( infile1 , 1 );

      //  2 --> means split mode, so load all into Xtrain
      load( infile2 , 2 );

      merge( outfile );
      
      return;
     
    }

  
  //
  // otherwise, always requires input data
  //
  
  const std::string infile = param.requires( "load" );


  //
  // split mode
  //

  if ( split_mode )
    {

      // 1 --> means split mode, so load all into Xtrain
      load( infile , 1 );
      
      std::string id_file1 = param.requires( "ids1" );
      std::string id_file2 = param.requires( "ids2" );
      std::string out1 = param.requires( "out1" );
      std::string out2 = param.requires( "out2" );
      if ( param.has( "vars" ) )
	{
	  std::set<std::string> v = param.strset( "vars" );
	  split( id_file1, id_file2, out1, out2 , &v );
	}
      else
	split( id_file1, id_file2, out1, out2 );
      return;
    }


  //
  // load observation types
  //
  
  attach_ids( param );
    
  //
  // load data 
  //

  load( infile );

  //
  // attach phenotypes
  //

  attach_phenotypes( param );

  //
  // prune training/validation datasets, if needed
  //

  prune();

  //
  // train
  //

  if ( train_mode )
    {
      train( param );

      save_model( param );
      
      return;
    }
  
  //
  // test 
  //

  if ( test_mode )
    {
      load_model( param );
      
      predict( param );
  
      SHAP( param );

      return;      
    }
  
  
  
  //
  // dump 
  //

  if ( dump_mode )
    {

      std::string ofile = Helper::expand( param.requires( "dump-training" ) );

      std::ofstream O1( ofile.c_str() , std::ios::out );

      logger << "  dumping training matrix to " << ofile << "\n";

      const int nrow = training_iids.size();
      const int ncol = vars.size();
      
      O1 << "IID\tID\tEID";
      for (int j=0; j<ncol; j++)
	O1 << "\t" << vars[j] ;
      O1 << "\n";
      for (int i=0; i<nrow; i++)
	{
	  O1 << training_iids[i] << "\t" << training_ids[i] << "\t" << training_eids[i];
	  for (int j=0; j<ncol; j++)
	    O1 << "\t" << Xtrain(i,j);
	  O1 << "\n";
	}

    }
  
  // expect:
  //   load -> binary data file
  //   model -> LGBM model file
  
  //   train   run as training mode (otherwise test)
  //     config
  //     iter
  //     train-ids
  //     test-ids

  //  (else, test mode)
  //     test-ids
  
  //  split
  

  
}

void massoc_t::load( const std::string & filename , const int force_destin )
{

  // special case: in test mode, if no pool of test IDs specified, assume it is everybody
  const bool all_test = mode == 2 && test_pool.size() == 0 ;

  if ( ! force_destin )
    {
      if ( ! all_test ) 
	if ( training_pool.size() == 0 && validation_pool.size() == 0 && test_pool.size() == 0 )
	  Helper::halt( "no training/validation/test obs specified... quitting" );
    }
  else
    {
      if ( training_pool.size() != 0 || validation_pool.size() != 0 || test_pool.size() != 0 )
        Helper::halt( "training/validation/test obs should not be specified in split/merge mode... quitting" );
    }
  
  // get counts of IDs
  int obs_train = 0 , obs_valid = 0 , obs_test = 0;
  
  // 0=ignore, 1=train, 2=valid, 3=test
  std::vector<int> home; 

  
  // load data  

  if ( ! Helper::fileExists( Helper::expand( filename ) ) )
    Helper::halt( "could not open " + filename );

  std::ifstream IN1( Helper::expand( filename ).c_str() , std::ios::binary | std::ios::in );
  
  // track # of IIDs actually used
  std::set<std::string> cnt_train, cnt_valid, cnt_test;
  
  // header (unless second file in a merge mode )
  if ( force_destin != 2 )
    vars.clear();

  // do not clear other IDs if in merge mode (i.e. load called twice)
  if ( force_destin == 0 )
    {
      training_iids.clear(); training_ids.clear(); training_eids.clear();
      validation_iids.clear(); validation_ids.clear(); validation_eids.clear();
      test_iids.clear(); test_ids.clear(); test_eids.clear();
    }
  else if ( force_destin == 1 )
    {
      training_iids.clear(); training_ids.clear(); training_eids.clear();
    }
  else if ( force_destin == 2 )
    {
      validation_iids.clear(); validation_ids.clear(); validation_eids.clear();
    }
  
  // we want to be able to read in concatenated files
  // so just get IDs on first pass
  
  while ( 1 )
    {  

      //
      // rows : observations in this set
      //
      
      const int nrow = Helper::bread_int( IN1 );

      // all done?
      if ( IN1.bad() || IN1.eof() ) break;
      
      
      for (int i=0; i<nrow; i++)
	{

	  // IID 
	  const std::string iid = Helper::bread_str( IN1 );

	  // class-level ID (e.g. spindle type)
	  const std::string id = Helper::bread_str( IN1 );

	  // event-level ID (e.g. spindle count, within type)
	  const std::string eid = Helper::bread_str( IN1 );

	  if ( force_destin == 1 || training_pool.find( iid ) != training_pool.end() )
	    {
	      home.push_back( 1 );
	      training_ids.push_back( id );
	      training_iids.push_back( iid );
	      training_eids.push_back( eid );
	      cnt_train.insert( iid );
	      ++obs_train;
	    }
	  else if ( force_destin == 2 || validation_pool.find( iid ) != validation_pool.end() )
	    {
	      home.push_back( 2 );
	      validation_ids.push_back( id );
	      validation_iids.push_back( iid );
	      validation_eids.push_back( eid );
	      cnt_valid.insert( iid );
	      ++obs_valid;
	    }
	  else if ( all_test || test_pool.find( iid ) != test_pool.end() )
	    {
	      home.push_back( 3 );
	      test_ids.push_back( id );
	      test_iids.push_back( iid );
	      test_eids.push_back( eid );
	      cnt_test.insert( iid );
	      ++obs_test;
	    }
	  else
	    {
	      home.push_back( 0 );
	    }
	  
	} // next obs


      //
      // now, variables
      //  if we've already read a subset, then need to check
      //  for alignment 
      //
      
      const int ncol = Helper::bread_int( IN1 );
      
      if ( vars.size() == 0 )
	{
	  for (int i=0; i<ncol; i++)
	    vars.push_back( Helper::bread_str( IN1 ) );
	}
      else
	{
	  if ( vars.size() != 0 && vars.size() != ncol )
	    Helper::halt( "mismatched # of columns in different subsets of " + filename );
	  
	  for (int i=0; i<ncol; i++)
            if ( vars[i] != Helper::bread_str( IN1 ) )
	      Helper::halt( "mismatched column header " + vars[i] );
	  
	}

      //
      // now skip the data
      //
      
      Helper::bskip_dbl( IN1 , nrow * ncol );

    }
  
  IN1.close();
  
  //
  // allocate storage
  //

  if ( force_destin == 1 )
    {
      Xtrain = Eigen::MatrixXd::Constant( obs_train, vars.size(),  NaN_value );
    }
  else if ( force_destin == 2 )
    {
      Xvalid = Eigen::MatrixXd::Constant( obs_valid, vars.size(),  NaN_value );
    }
  else
    {
      Xtrain = Eigen::MatrixXd::Constant( obs_train, vars.size(),  NaN_value );
      Xvalid = Eigen::MatrixXd::Constant( obs_valid, vars.size(),  NaN_value );
      Xtest  = Eigen::MatrixXd::Constant( obs_test , vars.size(),  NaN_value );
    }
  
  logger << "  reading feature matrices (" << vars.size() << " features in all cases):\n"
	 << "    - training   : " << obs_train << " obs from " << cnt_train.size() << " indivs\n"
  	 << "    - validation : " << obs_valid << " obs from " << cnt_valid.size() << " indivs\n"
	 << "    - test       : " << obs_test << " obs from " << cnt_test.size() << " indivs\n";
  

  //
  // re-read
  //

  std::ifstream IN2( Helper::expand( filename ).c_str() , std::ios::binary | std::ios::in );

  // reset counters
  obs_train = obs_valid = obs_test = 0;
  int obs_ignored = 0;

  int row = 0;
  
  while ( 1 )
    {  

      // rows - ignore now

      const int nrow = Helper::bread_int( IN2 );      

      if ( IN2.bad() || IN2.eof() ) break;
      
      for (int i=0; i<nrow; i++)
	{
	  const std::string iid = Helper::bread_str( IN2 );
	  const std::string id = Helper::bread_str( IN2 );
	  const std::string eid = Helper::bread_str( IN2 );
	}

      // cols - ignore now

      const int ncol = Helper::bread_int( IN2 );
      for (int i=0; i<ncol; i++)
	{
	  std::string dummy = Helper::bread_str( IN2 ) ;
	}

      // data
      
      for (int i=0; i<nrow; i++)
	{
	  if ( home[row] == 1 ) // training
	    {
	      for (int j=0; j<ncol; j++)
		Xtrain( obs_train , j ) = Helper::bread_dbl( IN2 ) ;
	      ++obs_train;
	    }
	  else if ( home[row] == 2 ) // validation
	    {
	      for (int j=0; j<ncol; j++)
		Xvalid( obs_valid , j ) = Helper::bread_dbl( IN2 ) ;
	      ++obs_valid;
	    }
	  else if ( home[row] == 3 ) // test
	    {
	      for (int j=0; j<ncol; j++)
		Xtest( obs_test , j ) = Helper::bread_dbl( IN2 ) ;
	      ++obs_test;
	    }
	  else // ignore
	    {
	      Helper::bskip_dbl( IN2 , ncol );
	      ++obs_ignored;
	    }
	  ++row;
	}
      
      // end of this subset... go back to see if any more blocks
    }
  
  // all done
  IN2.close();
  
  logger << " ... done\n";

  
  
}

void massoc_t::merge( const std::string & out_file )
{

  const int nrow = Xtrain.rows();
  if ( Xvalid.rows() != nrow ) Helper::halt( "load1 and load2 have different number of rows" );

  // datasets must correspond to the same individuals and events (i.e. IID and EID)
  //  that is - we allow the ID (spindle/channel class) to vary .. typically raw + band-pass filtered signal
  
  for (int i=0; i<nrow; i++)
    {
      if ( training_iids[i] != validation_iids[i] )
	Helper::halt( "load1 and load2 files do not have identical row structure" );
      if ( training_eids[i] != validation_eids[i] )
	Helper::halt( "load1 and load2 files do not have identical row structure" );      
    }
  
  // if here, okay to merge... 
  const int ncol = Xvalid.cols();
  const int ncol2 = 2 * ncol;
  
  // expand variable names
  std::vector<std::string> vars2 = vars;
  vars.resize( ncol2 );
  for (int j=ncol; j<ncol2; j++)
    vars[j] = vars2[ j-ncol ] + "_V2";

  // copy over data (into 'training')
  Xtrain.conservativeResize( Eigen::NoChange , ncol2 );
  
  for (int i=0; i<nrow; i++)
    {
      // ID
      training_ids[i] = training_ids[i] + "+" + validation_ids[i];
      // actual data
      for (int j=ncol; j<ncol2; j++) Xtrain(i,j) = Xvalid(i,j-ncol);
    }  

  // save the merged file

  logger << "  saving the merged dataset to " << out_file << "\n";
  
  save( training_iids, training_ids, training_eids, vars , Xtrain , out_file );
    
}


void massoc_t::split( const std::string & id_file1, const std::string & id_file2,
		      const std::string & out_file1, const std::string & out_file2 ,
		      const std::set<std::string> * xvars )
{

  // split always operates on 'Xtrain', arbitrarily
  //  i.e. split will load into Xtrain

  // training_iids[]  indiv-level IIDs
  // training_ids[]   obs-level IDs
  // vars[]
  // Xtrain[]
  
  // 0 = ignore , 1 = file1, 2 = file2 , 3 = both 
  std::vector<int> row( training_iids.size() );
  
  std::set<std::string> ids1, ids2;
  if ( ! Helper::fileExists( Helper::expand( id_file1 ) ) )
    Helper::halt( "could not open " + id_file1 );
  if ( ! Helper::fileExists( Helper::expand( id_file2 ) ) )
    Helper::halt( "could not open " + id_file2 );

  std::ifstream IN1( Helper::expand( id_file1 ).c_str() , std::ios::in );
  while ( 1 )
    {
      std::string l;
      IN1 >> l;
      if ( IN1.eof() || IN1.bad() ) break;
      if ( l == "" ) continue;      
      ids1.insert( l );
    }
  IN1.close();

  std::ifstream IN2( Helper::expand( id_file2 ).c_str() , std::ios::in );
  while	( 1 )
    {
      std::string l;
      IN2 >> l;
      if ( IN2.eof() || IN2.bad() ) break;
      if ( l == "" ) continue;      
      ids2.insert( l );
    }
  IN2.close();

  logger << "  read " << ids1.size() << " and " << ids2.size() << " IDs from "
	 << id_file1 << " and " << id_file2 << " respectively\n";

  std::vector<std::string> iids1, iids2;
  std::vector<std::string> rowids1, rowids2;
  std::vector<std::string> eids1, eids2;

  std::set<std::string> cnt_iids1, cnt_iids2, cnt_iids;
  
  int i1 = 0 , i2 = 0, i3 = 0;
  for (int i=0; i<training_iids.size(); i++)
    {
      cnt_iids.insert( training_iids[i] );

      // does the IID match?
      const bool m1 = ids1.find( training_iids[i] ) != ids1.end();
      const bool m2 = ids2.find( training_iids[i] ) != ids2.end();
      
      if ( m1 && m2 )
	{
	  row[i] = 3;
	  ++i1; ++i2; ++i3;
	  iids1.push_back(  training_iids[i] );
	  iids2.push_back(  training_iids[i] );
	  cnt_iids1.insert( training_iids[i] );
	  cnt_iids2.insert( training_iids[i] );
	  rowids1.push_back( training_ids[i] );
	  rowids2.push_back( training_ids[i] );
	  eids1.push_back( training_eids[i] );
	  eids2.push_back( training_eids[i] );
	}
      else if ( m1 )
	{
	  row[i] = 1;
	  ++i1;
	  iids1.push_back(  training_iids[i] );	  
	  cnt_iids1.insert( training_iids[i] );
	  rowids1.push_back( training_ids[i] );
	  eids1.push_back( training_eids[i] );
	}
      else if ( m2 )
	{
	  row[i] = 2;
	  ++i2;
	  iids2.push_back(  training_iids[i] );	  
	  cnt_iids2.insert( training_iids[i] );
	  rowids2.push_back( training_ids[i] );
	  eids2.push_back( training_eids[i] ); 
	}
      else
	row[i] = 0;
    }
  
  logger << "  original dataset: " << training_ids.size() << " observations on " << cnt_iids.size() << " individuals\n";
  logger << "        split 1 --> " << rowids1.size() << " observations on " << cnt_iids1.size()  << " individuals (" << out_file1 << ")\n";
  logger << "        split 2 --> " << rowids2.size() << " observations on " << cnt_iids2.size()  << " individuals (" << out_file2 << ")\n";
  if ( i3 )
    logger << "        (" << i3 << " observations will be written to both outputs)\n";
  
  //
  // split variables?
  //
  
  int ncol = vars.size();
  std::vector<bool> keep( ncol , true );
  int ncol_ret = ncol;
  std::vector<std::string> vars_ret = vars;
  if ( xvars != NULL )
    {
      vars_ret.clear();
      for (int j=0; j<ncol; j++)
	if ( xvars->find( vars[j] ) == xvars->end() )
	  {
	    keep[j] = false;
	    --ncol_ret;
	  }
	else
	  vars_ret.push_back( vars[j] );
    }
  
  if ( ncol_ret < ncol )
    logger << "  " << ncol_ret << " of " << ncol << " features will be extracted\n";
  
  
  //
  // write
  //
  
  Eigen::MatrixXd X1 = Eigen::MatrixXd::Zero( i1 , ncol_ret );
  Eigen::MatrixXd X2 = Eigen::MatrixXd::Zero( i2 , ncol_ret );
  
  i1 = 0; i2 = 0;
  
  for (int i=0; i<row.size(); i++)
    {
      if ( row[i] == 3 )
	{
	  X1.row( i1++ ) = Xtrain.row( i );
	  X2.row( i2++ ) = Xtrain.row( i );
	}
      else if ( row[i] == 2 )
	X2.row( i2++ ) = Xtrain.row( i );
      else if ( row[i] == 1 )
	X1.row( i1++ ) = Xtrain.row( i );
    }

  // write new files out
  save( iids1, rowids1 , eids1, vars_ret , X1 , out_file1 );
  save( iids2, rowids2 , eids2, vars_ret , X2 , out_file2 );
  
}

void massoc_t::save( const std::string & iid,
		     const std::vector<std::string> & rowids ,
		     const std::vector<std::string> & eids ,
		     const std::vector<std::string> & colids ,
		     const Eigen::MatrixXd & X , 
		     const std::string & filename )
{
  // expand single IID to vector
  std::vector<std::string> iids( rowids.size() , iid );
  save( iids , rowids , eids, colids , X , filename );
}


void massoc_t::save( const std::vector<std::string> & iids,
		     const std::vector<std::string> & rowids ,
		     const std::vector<std::string> & eids ,
		     const std::vector<std::string> & colids ,
		     const Eigen::MatrixXd & X , 
		     const std::string & filename )
{


  const int nrow = X.rows();
  const int ncol = X.cols();
  
  if ( colids.size() != ncol )
    Helper::halt( "mismatch in # of cols" );
  
  if ( rowids.size() != nrow )
    Helper::halt( "mismatch in # of rows" );

  if ( eids.size() != nrow )
    Helper::halt( "mismatch in # of rows" );
  
  if ( iids.size() != nrow )
    Helper::halt( "mismatch in # of rows" );
  
  // save data

  std::ofstream OUT1( Helper::expand( filename ).c_str() , std::ios::binary | std::ios::out );
  
  logger << " writing binary data matrix, " << ncol << " features, " << nrow << " observations\n";
  
  // rows 
  Helper::bwrite( OUT1, nrow ) ;
  for (int i=0; i<nrow; i++)
    {
      Helper::bwrite( OUT1, iids[i] ) ;  // IID
      Helper::bwrite( OUT1, rowids[i] ); // event class (e.g. CH.C3_F.11)
      Helper::bwrite( OUT1, eids[i] );   // event ID (1,2,3,...)
    }
  
  // cols
  Helper::bwrite( OUT1, ncol ) ;
  for (int i=0; i<ncol; i++)
    Helper::bwrite( OUT1 , colids[i] );
  
  // data
  for (int i=0; i<nrow; i++)
    for (int j=0; j<ncol; j++)
      Helper::bwrite( OUT1 , X(i,j) );      
  
  OUT1.close();
  
}



// when called internally (e.g. from TLOCK, where we just have one big feature matrix)

massoc_t::massoc_t( const std::string & iid,
		    const std::vector<std::string> & rowids ,
		    const std::vector<std::string> & eids ,
		    const std::vector<std::string> & colids ,
		    const Data::Matrix<double> & XX ,
		    const std::string & filename )
{

  // yes, dumb copy, I know
  // also note::: whilst coming from TLOCK, XX needs to be transposed --> X

  const int nrow = XX.dim2();
  const int ncol = XX.dim1();
  
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero( nrow, ncol );
  for (int i=0; i<nrow; i++)
    for (int j=0; j<ncol; j++)
      X(i,j) = XX(j,i); // nb. transpose
  
  save( iid, rowids, eids, colids, X , filename );
  
}

  
// load phenotypes/labels
void massoc_t::attach_phenotypes( param_t & param )
{
  
  Ytrain.resize( Xtrain.rows() , NaN_value );
  Yvalid.resize( Xvalid.rows() , NaN_value );
  Ytest.resize( Xtest.rows() , NaN_value );
  
  // by default match on full IID+ID+EID (i.e. this is expected in the vars file)
  // however, if iid-vars=T, then match on IID only (from ID in vars).  
  //   i.e. here we apply individual (not event) level labels
  
  const bool iid_match = param.yesno( "iid-vars" );
  
  // known phenotypes optional for test dataset                                                                                                 
  if ( ! param.has( "phe" ) ) return;
  
  phenotype_label = param.value( "phe" );

  double y;

  int obs_train = 0 , obs_valid = 0 , obs_test = 0;
  // training data
  for (int i=0; i<training_ids.size(); i++)
    {
      const std::string id = iid_match ? training_iids[i] : training_iids[i] + "_" + training_ids[i] + "_" + training_eids[i];	
      if ( cmd_t::pull_ivar( id , phenotype_label , &y ) ) { Ytrain[i] = y; ++obs_train; }
      
    }
  
  // validation data
  for (int i=0; i<validation_ids.size(); i++)
    {
      const std::string id = iid_match ? validation_iids[i] : validation_iids[i] + "_" + validation_ids[i] + "_" + validation_eids[i];
      if ( cmd_t::pull_ivar( id , phenotype_label , &y ) ) { Yvalid[i] = y; ++obs_valid; }      
    }

  
  // test data
  for (int i=0; i<test_ids.size(); i++)
    {
      const std::string id = iid_match ? test_iids[i] : test_iids[i] + "_" + test_ids[i] + "_" + test_eids[i];
      if ( cmd_t::pull_ivar( id , phenotype_label , &y ) ) { Ytest[i] = y; ++obs_test; }       
    }  
  
  logger << "  attached " << phenotype_label << " for " 
	 << obs_train << " (of " << training_ids.size()  << ") training, "
	 << obs_valid << " (of " << validation_ids.size() << ") validation, and "
	 << obs_test  << " (of " << test_ids.size() << ") test observation\n";
  
}

// prune trainers/validation data w/ a missing label
void massoc_t::prune()
{
 
  // trainers
  
  const int ni_train = Ytrain.size() ;
  std::vector<bool> miss_train( ni_train );
  int obs_train = 0;
  for (int i=0; i<ni_train; i++)
    {
      miss_train[i] = std::isnan( Ytrain[i] ) ;
      if ( ! miss_train[i] ) ++obs_train;
    }
  
  if ( obs_train < ni_train )
    {
      logger << "  pruning train dataset from " << ni_train << " to " << obs_train << " based on missing/NA labels\n";
      prune1( obs_train , miss_train , &training_iids, &training_ids, &training_eids, &Xtrain, &Ytrain );
    }

  
  // validation
  
  const int ni_valid = Yvalid.size() ;
  std::vector<bool> miss_valid( ni_valid );
  int obs_valid = 0;
  for (int i=0; i<ni_valid; i++)
    {
      miss_valid[i] = std::isnan( Yvalid[i] ) ;
      if ( ! miss_valid[i] ) ++obs_valid;
    }
  
  if ( obs_valid < ni_valid )
    {
      logger << "  pruning validation dataset from " << ni_valid << " to " << obs_valid << " based on missing/NA labels\n";
      prune1( obs_valid , miss_valid , &validation_iids, &validation_ids, &validation_eids, &Xvalid, &Yvalid );
    }


}


void massoc_t::prune1( const int n, const std::vector<bool> & missing , 
		       std::vector<std::string> * iids, 
		       std::vector<std::string> * ids, 
		       std::vector<std::string> * eids, 
		       Eigen::MatrixXd * X, 
		       std::vector<double> * Y )
{

  const int n0 = iids->size();

  if ( ids->size() != n0 || 
       eids->size() != n0  ||
       X->rows() != n0 ||
       Y->size() != n0 ||
       missing.size() != n0 )
    Helper::halt( "internal error in prune()" );
  
  // allocate copy storage
  std::vector<std::string> iids2( n );
  std::vector<std::string> ids2( n );
  std::vector<std::string> eids2( n );
  Eigen::MatrixXd X2 = Eigen::MatrixXd::Zero( n , X->cols() );
  std::vector<double> Y2( n );
  
  int c = 0;
  for (int i=0; i<n0; i++)
    {
      if ( ! missing[i] ) 
	{
	  iids2[c] = (*iids)[i];
	  ids2[c] = (*ids)[i];
	  eids2[c] = (*eids)[i];
	  X2.row(c) = X->row(i);
	  Y2[c] = (*Y)[i];
	  ++c;
	}
    }

  // copy back
  *iids = iids2;
  *ids = ids2;
  *eids = eids2;
  *X = X2;
  *Y = Y2;
}




// load train/valid/test status
void massoc_t::attach_ids( param_t & param )
{

  training_pool.clear();
  validation_pool.clear();
  test_pool.clear();

  //
  // Training IDs
  //

  if ( param.has( "train-ids" ) )
    {
      const std::string filename = Helper::expand( param.value( "train-ids" ) ) ;
      
      if ( ! Helper::fileExists( filename ) )
	Helper::halt( "could not open " + filename );
      
      std::ifstream IN1( filename.c_str(), std::ios::in );
      while ( 1 )
	{
	  std::string l;
	  IN1 >> l;
	  if ( IN1.eof() || IN1.bad() ) break;
	  if ( l == "" ) continue;
	  training_pool.insert( l );
	}
      IN1.close();
    }
  

  //
  // Validation IDs
  //

  if ( param.has( "valid-ids" ) )
    {
      const std::string filename = Helper::expand( param.value( "valid-ids" ) ) ;
      
      if ( ! Helper::fileExists( filename ) )
	Helper::halt( "could not open " + filename );
      
      std::ifstream IN1( filename.c_str(), std::ios::in );
      while ( 1 )
	{
	  std::string l;
	  IN1 >> l;	  
	  if ( IN1.eof() || IN1.bad() ) break;
	  if ( l == "" ) continue;
	  validation_pool.insert( l );
	}
      IN1.close();
    }


  //
  // Test IDs
  //
  
  if ( param.has( "test-ids" ) )
    {
      const std::string filename = Helper::expand( param.value( "test-ids" ) ) ;
      
      if ( ! Helper::fileExists( filename ) )
        Helper::halt( "could not open " + filename );
      
      std::ifstream IN1( filename.c_str(), std::ios::in );
      while ( 1 )
        {
          std::string l;
          IN1 >> l;
          if ( IN1.eof() || IN1.bad() ) break;
	  if ( l == "" ) continue;          
          test_pool.insert( l );
        }
      IN1.close();
    }

  logger << "  read " << training_pool.size() << " training IDs, "
	 << validation_pool.size() << " validation IDs, and "
	 << test_pool.size() << " test IDs\n";
  
}


//
// train model
//

void massoc_t::train( param_t & param )
{
  // configuration required  
  lgbm.load_config( param.requires( "config" ) );

  // std::cout << " Y train -  \n";
  // for (int i=0; i<Ytrain.size(); i++) std::cout << Ytrain[i] << "\n";
  
  // attach data
  lgbm.attach_training_matrix( Xtrain );
  lgbm.attach_training_qts( Ytrain );

  if ( validation_ids.size() > 0 )
    {
      lgbm.attach_validation_matrix( Xvalid );
      lgbm.attach_validation_qts( Yvalid );
    }
  
  // (TODO) attach weights?
  // ...

  // specify the number of iterations? (default = 100) 
  lgbm.n_iterations = param.has( "iter" ) ? param.requires_int( "iter" ) : 100 ; 

  // train model
  lgbm.create_booster( true );

}

void massoc_t::save_model( param_t & param )
{
  //
  // write out trained LGBM model
  //

  const std::string model_file = param.requires( "model" );
  
  lgbm.save_model( model_file );

}


//
// prediction
//

void massoc_t::load_model( param_t & param )
{
  const std::string model_file = param.requires( "model" );
  
  lgbm.load_model( model_file );
  
  logger << "  read LGBM model file from " << model_file << "\n";

}
  
void massoc_t::predict( param_t & param )
{

  // only go up to iteration 'iter'?  ( 0 implies all )
  const int iter = param.has( "iter" ) ? param.requires_int( "iter" ) : 0 ;
  
  Eigen::MatrixXd Y = lgbm.predict( Xtest , iter );

  const int n = Y.rows();

  if ( test_ids.size() != n )
    Helper::halt( "internal error in predict()" );

  writer.id( "." , "." );

  bool allobs = true;
  
  for (int i=0;i<n;i++)
    {
      writer.id( test_iids[i] + "_" + test_ids[i] + "_" + test_eids[i] , "." );

	   
      writer.value( "IID" , test_iids[i] ); // indiv
      writer.value( "TID" , test_ids[i] );  // type

      int eid = 0;
      if ( Helper::str2int( test_eids[i] , &eid ) )
	writer.value( "EID" , eid ); // should always work
      else
	writer.value( "EID" , test_eids[i] ); // but in case not...

      if ( Ytest.size() != 0 && Helper::realnum( Ytest[i] ) )
	writer.value( "OBS" , Ytest[i] );
      else
	allobs = false;
      
      writer.value( "PRD" , Y(i,0) );
    }
    
  
  writer.id( "." , "." );

  //
  // get a correlation between predicted/observed (overall) 
  //

  if ( allobs )
    {
      // change these...
      double r = Statistics::correlation( Ytest, eigen_ops::copy_array( Y.col(0) ) );
      writer.value( "R" , r );
    }


}



void massoc_t::SHAP( param_t & param )
{

  // only go up to iteration 'iter'?  ( 0 implies all )
  const int iter = param.has( "iter" ) ? param.requires_int( "iter" ) : 0 ;
  
  Eigen::MatrixXd S = lgbm.SHAP_values( Xtest , iter );
  
  const int n = S.rows();
  const int nv = S.cols() - 1; // nb last col is expected value

  if ( test_ids.size() != n )
    Helper::halt( "internal error in predict()" );

  //  std::cout << " nv = " << nv << " " << varlist.size() << "\n";
    
  if ( nv != vars.size() )
    Helper::halt( "internal error in predict(), varlist size" );
  
  writer.id( "." , "." );
  
  Eigen::VectorXd M = S.cwiseAbs().colwise().mean();

  if ( M.size() != nv + 1 )
    Helper::halt( "internal error in SHAP" );

  for (int j=0; j<nv; j++)
    {
      writer.level( vars[j] , "VAR" );
      writer.value( "SHAP" , M[j] );
    }
  writer.unlevel( "VAR" );

  // indiv level output
  
  if ( param.has( "verbose" ) ) 
    {
      for (int i=0;i<n;i++)
	{
	  writer.id( test_iids[i] + "_" + test_ids[i] + "_" + test_eids[i] , "." );
	  
	  for (int j=0; j<nv; j++)
	    {
	      writer.level( vars[j] , "VAR" );
	      writer.value( "SHAP" , S(i,j) );
	    }
	  writer.unlevel( "VAR" );
	}
      writer.id( "." , "." );
    }
 

}


#endif
