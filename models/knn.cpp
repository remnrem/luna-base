
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
  
void model_knn_t::load( const std::string & f , int * rows , int * cols )
{
  std::cout << "X.rows() = " << X.rows() << "\n";
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
      X(i,j) = x;
      ++j;
      if ( j == ncol )
	{
	  ++i;
	  j = 0;
	}
    }

  //  std::cout << "X\n" << X << "\n";
  
}
  
  
Eigen::VectorXd model_knn_t::impute( const Eigen::VectorXd & f ,
				     const std::vector<bool> & missing ,
				     int k )
{

  if ( X.rows() == 0 ) return f;

  // ** assumes that reference is complete - no missing )

  const int n = f.size();

  std::vector<int> slots;
  for (int i=0; i<n; i++)
    if ( ! missing[i] ) slots.push_back( i );

  const int n1 = slots.size();
  const int ni = X.rows();
  
  Eigen::MatrixXd X1 = Eigen::MatrixXd::Zero( ni , n1 );
  for (int i=0; i<slots.size(); i++)
    X1.col(i) = X.col(slots[i]);

  std::map<double,int> nearest;
  for (int i=0; i<ni; i++)
    nearest[ ( f - X1.row(i) ).norm() ] = i; 

  //
  std::vector<double> imps( n1 , 0 );

  int j = 0;
  std::map<double,int>::const_iterator nn = nearest.begin();
  while ( nn != nearest.end() )
    {
      if ( j == k ) break;
      int i = nn->second;
      for (int q=0; q<n1; q++)
	imps[q] += X(i,slots[q]);
      ++j;
    }

  // put means back in
  Eigen::VectorXd f2 = f;  
  for (int q=0; q<n1; q++)
    f2[slots[q]] = imps[q] / (double)k;
  
  // all done
  return f2;
  
}
  
 
