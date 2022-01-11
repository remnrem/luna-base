

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

#include "stats/qda.h"
#include "stats/statistics.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "miscmath/miscmath.h"
#include "stats/eigen_ops.h"

#include <cmath>
#include <iostream>
#include <fstream>

extern logger_t logger;


qda_model_t qda_t::fit( const bool flat_priors )
{

  qda_model_t model;

  //
  // Some observatons can be 'missing' and will be skipped 
  // in generating the model
  //
  
  int nm = 0;
  std::vector<bool> mi( X.rows() , false );
  for (int i=0; i<y.size(); i++)
    if ( y[i] == missing ) 
      { ++nm; mi[i] = true; }
  
  const int n = X.rows() - nm;
  const int p = X.cols();
  
  if ( n < 3 ) 
    {
      logger << "  *** not enough nonmissing obs for QDA\n" ;
      model.valid = false;
      return model;      
    }

  if ( nm )
    {      
      std::vector<std::string> y2 = y;
      Eigen::MatrixXd X2 = X;
      
      y.resize( n );

      X.resize( n , p );
      int n0 = y2.size();
      int c = 0;
      for (int i=0; i<n0; i++)
	{
	  if ( mi[i] ) continue;
	  y[c] = y2[i];
	  for (int j=0; j<p; j++)
	    X(c,j) = X2(i,j);
	  ++c;
	}
    }

  //
  // At this point, any missing values will have been removed, and we can 
  // fit the QDA model w/out worrying about missing values
  //


  //
  // count group label
  //

  std::map<std::string,int> counts;
  for (int i=0;i<y.size();i++) 
    counts[ y[i] ]++;

  
  //
  // assign index to groups (integer 0, 1, ...)
  //
  
  std::map<std::string,int> gidx;
  std::map<std::string,int>::const_iterator cc = counts.begin();
  while ( cc != counts.end() )
    {      
      const int sz = gidx.size();
      gidx[ cc->first ] = sz;
      ++cc;
    }
  
  std::vector<int> yi( n );
  for (int i=0;i<n;i++) 
    yi[i] = gidx[y[i]];    
    
  //
  // number of classes
  //

  const int ng = counts.size();
  
  if ( ng == 1 ) 
    {
      logger << "  *** no variation in group labels for QDA\n";
      model.valid = false;
      return model;      
    }

  //
  // priors
  //
  
  ArrayXd prior( ng );
  ArrayXd nc( ng );
  
  int cidx = 0;
  cc = counts.begin();
  while ( cc != counts.end() )
    {
    
      if ( cc->second < p + 1 ) 
	{
	  logger << "  ** group size too small for QDA\n";
	  model.valid = false;
	  return model;
	}

      if ( flat_priors )
	prior[ cidx ] = 1.0 / (double) counts.size() ;
      else
	prior[ cidx ] = cc->second / (double)n ;

      // store count N
      nc[ cidx ] = cc->second;
      
      ++cc; 
      ++cidx;
    }
  

  //
  // group means () ng x p ( groups x predictors )
  //
  
  Eigen::MatrixXd group_means = Eigen::MatrixXd::Zero( ng , p );
  for ( int i = 0 ; i < n ; i++ )
    for ( int j = 0 ; j < p ; j++ )
      group_means( yi[i], j ) += X(i,j);
  
  for ( int i = 0 ; i < ng ; i++ )
    for ( int j = 0 ; j < p ; j++ )
      group_means( i, j ) /= n * prior[i];
  
  
  //
  // Perform each class at a time
  //

  std::vector<Eigen::MatrixXd> scaling( ng );
  std::vector<double> ldet( ng );
  
  for (int i = 0 ; i < ng ; i++ )
    {      
      
      const double sqrt_nk = sqrt( nc[i]  - 1 );
      
      // mean-centered, group-specific X matrix ('X1')

      Eigen::MatrixXd X1 = Eigen::MatrixXd::Zero( nc[i] , p );
      int idx = 0;
      for (int j = 0 ; j < n ; j++)
	if ( yi[j] == i )
	  {
	    X1.row(idx) = ( X.row(j) - group_means.row( i ) ) / sqrt_nk ; 
	    ++idx;
	  }

      //
      // check within-group variability
      //
      
      for (int j=0; j<p; j++)
	{
	  double s = eigen_ops::sdev( X1.col(j) );
	  
	  if ( s < tol )
	    {
	      model.valid = false;
	      model.errmsg = "variable " + Helper::int2str(j) + " is constant within group ";
	      return model;
	    }      
	}
      
      
      //
      // QR decomposition to get scaling matrix (to whiten X)
      //

      Eigen::HouseholderQR<Eigen::MatrixXd> qr( X1 );      
      
      Eigen::MatrixXd R = qr.matrixQR().topRows( p );

      scaling[i] = R.triangularView<Eigen::Upper>().solve( Eigen::MatrixXd::Identity( p , p ) );
      
      // log abs det 
      ldet[i] = R.diagonal().array().abs().log().sum() * 2;
      	        
    }

  //
  // Report
  //

  model.valid = true;
  model.prior = prior;
  model.counts = counts;
  model.rows = nc;
  model.means = group_means;
  model.scaling = scaling;
  model.ldet = ldet;
  model.n = n;
  model.labels.clear();
  
  std::map<std::string,int>::const_iterator ll = counts.begin();
  while ( ll != counts.end() )
    {
      model.labels.push_back( ll->first );
      ++ll;
    }

  return model;
}
  

qda_posteriors_t qda_t::predict( const qda_model_t & model , const Eigen::MatrixXd & X )
{
  
  if ( ! model.valid ) 
    Helper::halt( "internal error: QDA predict() is being passed an invalid model" );

  const int p = X.cols();
  const int n = X.rows();

  //  std::cout << "qda_t::predict() " << p << " " << n << "\n";
  
  if ( p != model.means.cols() )
    Helper::halt( "wrong number of columns in qda_t::predict(): expecting "
		  + Helper::int2str((int) model.means.cols() ) + " but found " + Helper::int2str( p ) );  
  
  const int ng = model.prior.size();

  // std::cout << "ng = " << ng << "\n";

  // std::cout << " means " << model.means << "\n\n"
  // 	    << " ldet " << model.ldet[0] << " " << model.ldet[1] << " .. " << model.ldet.size() << "\n";
  //
  // Predict 
  //

  Eigen::MatrixXd D = Eigen::MatrixXd::Zero( n , ng );
  
  for (int g=0; g<ng; g++)
    {
      // dev        <- ((x - matrix(object$means[i,  ], nrow(x), ncol(x), byrow = TRUE)) %*% object$scaling[,,i])
      Eigen::MatrixXd X1 = X.rowwise() - model.means.row(g) ;
      X1 = X1 * model.scaling[g] ; 

      // // dist[, i]  <- 0.5 * rowSums(dev^2) + 0.5 * object$ldet[i] - log(prior[i])
      D.col(g) = 0.5 * X1.array().square().rowwise().sum() + 0.5 * model.ldet[g] - log(model.prior[g]) ;
             
    }

  // normalize, 
  //dist <- exp( -(dist - apply(dist, 1L, min, na.rm=TRUE)))
  Eigen::VectorXd m = D.rowwise().minCoeff();
  D = D.colwise() - m ;  
  D = -1 * D;
  D = D.array().exp() ;

  // as probs.
  // posterior <- dist/drop(dist %*% rep(1, ngroup))
  Eigen::VectorXd s = D.rowwise().sum() ;
  
  for (int g=0;g<ng; g++)
    D.col(g) = D.col(g).cwiseQuotient( s );


  //
  // labels
  //

  std::vector<int> cli( n );
  std::vector<std::string> cl( n );
  for (int i=0;i<n;i++)
    {
      int mxi = 0;
      for (int k=1;k<ng;k++) if ( D(i,k) > D(i,mxi) ) mxi = k;
      cli[i] = mxi;
      cl[i] = model.labels[ mxi ];
      //std::cout << " cl " << cl[i] << "\n";
    }
  
  //
  // Construct final return object
  //
  			      
  qda_posteriors_t pp;

  // posterior probs
  pp.pp = D;
  
  // most likely class assignment (idx)
  pp.cli = cli;
  
  // as above, but using string label from model
  pp.cl = cl;
  
  return pp;

} 


void qda_model_t::write( const std::string & filename ) const
{
  if ( ! valid ) Helper::halt( "cannot write an invalid model" );
  std::ofstream O1( Helper::expand( filename ).c_str() , std::ios::out );
  
  O1 << "QDA\n";
  O1 << "ng: " << prior.size() << "\n";
  O1 << "nf: " << means.cols() << "\n";
  
  O1 << "priors:";
  for (int i=0;i<prior.size();i++) O1 << " " << prior[i] ;
  O1 << "\n";

  O1 << "rows:";
  for (int i=0;i<rows.size();i++) O1 << " " << rows[i] ;
  O1 << "\n";

  O1 << "counts:";
  std::map<std::string,int>::const_iterator ii = counts.begin();
  while ( ii != counts.end() ) 
    {
      O1 << " " << ii->first << " " << ii->second;
      ++ii;
    }
  O1 << "\n";

  O1 << "means:\n" << means
     << "\n";
  
  O1 << "scaling:\n";
  for (int i=0; i<scaling.size(); i++)
    O1 << scaling[i] << "\n";
  
  O1 << "ldet:";
  for (int i=0;i<ldet.size();i++) O1 << " " << ldet[i] ;
  O1 << "\n";

  O1 << "n: " << n << "\n";
  
  O1 << "labels:";
  for (int i=0;i<labels.size();i++) O1 << " " << labels[i] ;
  O1 << "\n";
  
  O1.close();
}


void qda_model_t::read( const std::string & filename )
{
  if ( ! Helper::fileExists( Helper::expand( filename ) ) )
    Helper::halt( "could not open " + filename );

  valid = true;
  errmsg = "";
  
  std::ifstream I1( Helper::expand( filename ).c_str() , std::ios::in );
  std::string dummy;

  int ng, nf;
  
  I1 >> dummy; // QDA
  
  I1 >> dummy >> ng 
     >> dummy >> nf;
  
  // size up 
  prior.resize( ng );
  rows.resize( ng );
  means.resize( ng , nf );
  scaling.resize( ng );
  for (int i=0;i<ng;i++) scaling[i].resize( nf, nf );
  ldet.resize( ng );
  labels.resize( ng );
  
  // read in 
  I1 >> dummy;
  for (int i=0;i<prior.size();i++) I1 >> prior[i] ;

  I1 >> dummy;
  for (int i=0;i<rows.size();i++) I1 >> rows[i] ;
  
  I1 >> dummy;
  counts.clear();
  for (int i=0;i<ng;i++)
    {
      std::string l;
      int c;
      I1 >> l >> c;
      counts[l] = c ;
    }

  I1 >> dummy;
  for (int i=0;i<ng;i++) 
    for (int j=0;j<nf;j++)
      I1 >> means(i,j);

  I1 >> dummy;
  for (int i=0;i<ng;i++) 
    for (int j=0;j<nf;j++)
      for (int k=0;k<nf;k++)
	I1 >> scaling[i](j,k);

  I1 >> dummy;
  for (int i=0;i<ldet.size();i++) I1 >> ldet[i] ;

  I1 >> dummy;
  I1 >> n;
  
  I1 >> dummy;
  for (int i=0;i<labels.size();i++) I1 >> labels[i] ;
    
  I1.close();


}
