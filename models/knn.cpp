
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

#include "models/knn.h"
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include "helper/helper.h"
#include "helper/logger.h"

extern logger_t logger;

void model_knn_t::clear()
{
  X = Eigen::MatrixXd::Zero(0,0);
}
  
void model_knn_t::load( const std::string & f ,
			const std::vector<std::string> & header , 
			int * rows , int * cols )
{

  // only load once
  if ( X.rows() != 0 )
    {
      if ( rows != NULL ) 
	*rows = X.rows();
      if ( cols != NULL ) 
	*cols = X.cols();
      return;
    }

  const std::string & filename = Helper::expand( f );

  if ( ! Helper::fileExists( filename ) )
    Helper::halt( "could not open " + filename );

  int nrow = 0 , ncol = 0;

  std::ifstream IN1( filename.c_str() , std::ios::in );
  std::string comm;

  // first row:
  // # nrows ncols
  IN1 >> comm >> nrow >> ncol;

  logger << "  creating " << nrow << " x " << ncol << " reference feature matrix from "
	 << filename << "\n";

  // second row: # column labels
  std::vector<std::string> labels( ncol );
  IN1 >> comm;  // another #
  for (int i=0; i<ncol; i++)
    IN1 >> labels[i];

  // check/get slot alignment
  if ( header.size() != ncol )
    Helper::halt( "expecting " + Helper::int2str( (int)header.size() ) + " columns but found " + Helper::int2str( ncol ) + " in " + filename );
  
  std::map<std::string,int> headerset;
  for (int i=0; i<ncol; i++)
    headerset[ header[i] ] = i;
  
  bool okay = true;
  bool order = true;
  for (int i=0; i<ncol; i++)
    {
      if ( headerset.find( labels[i] ) == headerset.end() )
	{
	  logger << "  could not find " << labels[i]  << " in the model specifcation\n";
	  okay = false;
	}
      if ( labels[i] != header[i] ) order = false;
    }
  
  if ( ! okay )    
    Helper::halt( "mismatch of column label count between model and data files " + filename );

  std::vector<int> slot(ncol);
  if ( order )
    for (int i=0; i<ncol; i++) slot[i] = i;
  else
    for (int i=0; i<ncol; i++) slot[i] = headerset[ labels[i] ] ;
    
  // now the data...
  X = Eigen::MatrixXd::Zero( nrow , ncol );
  
  if ( rows != NULL ) 
    *rows = X.rows();
  if ( cols != NULL ) 
    *cols = X.cols();

  int i = 0 , j = 0;
  while ( 1 )
    {
      double x ;
      
      IN1 >> x;
      if ( IN1.eof() ) break;
      
      // insert into the correct slot
      X(i,slot[j]) = x;
      
      ++j;
      if ( j == ncol )
	{
	  ++i;
	  j = 0;
	}
    }

  //    std::cout << "X\n" << X << "\n";
  
}
  
  
Eigen::VectorXd model_knn_t::impute( const Eigen::VectorXd & f ,
				     const std::vector<bool> & missing )

{

  if ( X.rows() == 0 ) return f;
  
  // ** assumes that reference is complete - no missing )
  
  const int n = f.size();
  
  // slots1 : has data - use to align
  // slots0 : missing data - need to fill in

  std::vector<int> slots1, slots0;
  
  for (int i=0; i<n; i++)
    if ( ! missing[i] )
      slots1.push_back( i );
    else
      slots0.push_back( i );
  
  const int n1 = slots1.size();
  const int n0 = slots0.size();
  const int ni = X.rows();

  // subset
  Eigen::MatrixXd X1 = Eigen::MatrixXd::Zero( ni , n1 );
  Eigen::VectorXd f1 = Eigen::VectorXd::Zero( n1 );
  for (int i=0; i<n1; i++)
    {
      X1.col(i) = X.col(slots1[i]);
      f1[i] = f[slots1[i]];
    }

  // rank by Euclidean distance
  std::map<double,int> nearest;
  for (int i=0; i<ni; i++)
    nearest[ ( f1 - X1.row(i).transpose() ).norm() ] = i; 

  // compile full vectors
  std::vector<double> imps( n , 0 );
  
  int j = 0;
  std::map<double,int>::const_iterator nn = nearest.begin();
  while ( nn != nearest.end() )
    {
      if ( j == k ) break;
      int i = nn->second;
      //      std::cout << " person " << j << " id = " << i << ", d = " << nn->first << "\n";
      for (int q=0; q<n; q++)
	imps[q] += X(i,q);
      ++j;
      ++nn;
    }

  // for (int q=0; q<n; q++)
  //   std::cout << " q/n = " << q << " " << imps[q] << "\n";
  // for (int q=0; q<n0; q++)
  //   std::cout << " q/n0 = " << q << " " << slots0[q] << " " << imps[slots0[q]] << "\n";

  
  // update feature vector for slots0 holes
  Eigen::VectorXd f2 = f;  
  for (int q=0; q<n0; q++)
    f2[slots0[q]] = imps[slots0[q]] / (double)k;
    
  // all done
  return f2;
  
}
  
 

Eigen::VectorXd model_knn_t::distance( const Eigen::VectorXd & f )				       
{
  
  // drop each term, one at a time and check
  // (only if non-missing)
  // no missing data allowed
  
  const int nt = f.size();
  const int ni = X.rows();
  
  if ( nt != X.cols() )
    Helper::halt( "feature vector does not align with training data" );

  // get means
  Eigen::VectorXd means = X.colwise().mean();
  
  // get SDs
  Eigen::Array<double, 1, Eigen::Dynamic> sds = ((X.rowwise() - X.colwise().mean()).array().square().colwise().sum()/(ni-1)).sqrt();

  Eigen::VectorXd V = Eigen::VectorXd::Zero( nt );
    
  // drop each feature, one at a time
  for (int i=0; i<nt; i++)
    {      
      std::vector<bool> missing( nt , false );
      missing[i] = true;

      // drop value, then re-impute
      Eigen::VectorXd I = impute( f , missing );

      // as we require inputs to be standardized, this gives how many SD units from the expected
      Eigen::VectorXd D = f - I;
      
      V[i] = D[i];      
    }
  
  return V;
  
}
