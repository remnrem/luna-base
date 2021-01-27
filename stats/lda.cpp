
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

#include "stats/lda.h"
#include "stats/statistics.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "miscmath/miscmath.h"

#include <cmath>
#include <iostream>

extern logger_t logger;


lda_model_t lda_t::fit( const bool flat_priors )
{

  lda_model_t model;

  const int n = X.rows();
  const int p = X.cols();

  //
  // count group labels
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
  for (int i=0;i<n;i++) { 
    yi[i] = gidx[y[i]];    
  }

  //
  // number of classes
  //

  const int ng = counts.size();
  
  if ( ng == 1 ) 
    Helper::halt( "no variation in group labels in lda_t::fit()" );

  //
  // priors
  //

  ArrayXd prior( ng );

  int cidx = 0;
  cc = counts.begin();
  while ( cc != counts.end() )
    {
      if ( flat_priors )
	prior[ cidx ] = 1.0 / (double) counts.size() ;
      else
	prior[ cidx ] = cc->second / (double)n ;
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
  // adjust X by group mean; get variance of each measure
  // (i.e. looking for within-group variability)
  //

  std::vector<double> f1(p);
  
  for (int j=0;j<p;j++)
    {
      std::vector<double> xa(n);
      for (int i=0;i<n;i++) 
	xa[i] = X(i,j) - group_means( yi[i] , j ) ;
      
      f1[j] = MiscMath::sdev( xa );
      
      if ( f1[j] < tol ) {
	model.valid = false;
	model.errmsg = "variable " + Helper::int2str(j) + " is constant within group ";
	return model;
      }      
    }

  
  // scaling matrix (diagonal)... 
  Eigen::MatrixXd scaling = Eigen::MatrixXd::Zero( p , p );
  for (int j=0;j<p;j++) 
    scaling(j,j) = 1.0 / f1[j] ;
  
  double fac = 1.0 / double( n - ng );
  
  double sqrt_fac = sqrt( fac );

  //
  // scaled X
  //

  Eigen::MatrixXd X1( n , p );
  for (int i=0;i<n;i++)
    for (int j=0;j<p;j++)      
      X1(i,j) = sqrt_fac * ( X(i,j) - group_means( yi[i] , j ) ) * scaling(j,j); 
  
  //
  // SVD
  //

  // bool okay = Statistics::svdcmp( X1 , W , V );
  // if ( ! okay ) Helper::halt( "problem in lda_t, initial SVD\n" );

  Eigen::BDCSVD<Eigen::MatrixXd> svd( X1 , Eigen::ComputeThinU | Eigen::ComputeThinV );  
  Eigen::MatrixXd V = svd.matrixV();
  ArrayXd W = svd.singularValues();
  
  // components are already sorted in decreasing order
  
  // X.s <- svd(X1, nu = 0L)
  
  //
  // rank
  //
  
  int rank = 0;
  for (int j=0;j<p;j++) 
    if ( W(j) > tol ) ++rank;
  
  if ( rank == 0 ) {
    model.valid = false;
    model.errmsg = "problem with collinearity/constant values in input data" ;
    return model;
  }

  if ( rank < p ) 
    logger << " warning... rank < p\n"; 


  //
  // get new scaling matrix p x rank , possibly of rank < p
  // scaling <- scaling %*% X.s$v[, 1L:rank] %*% diag(1/X.s$d[1L:rank],,rank)
  //

  Eigen::MatrixXd scaling2 = Eigen::MatrixXd::Zero( p , rank );
  for (int i=0;i<p;i++)
    for (int j=0;j<rank;j++)
      for (int k=0;k<p;k++)
	scaling2(i,j) += scaling(i,k) * V(k,j);

  for (int i=0;i<p;i++)
    for (int j=0;j<rank;j++)
      scaling2(i,j) /= W[j];

  //
  // mean values of each measure
  //
  
  ArrayXd xbar( p );
  for (int j=0;j<p;j++) xbar[j] = 0;

  for (int i=0;i<ng;i++)
    for (int j=0;j<p;j++)
      xbar[j] += prior[i] * group_means(i,j);
  
  
  fac = 1.0/(double)(ng - 1);
  
  
  //X <- sqrt((n * prior)*fac) * scale(group.means, center = xbar, scale = FALSE) %*% scaling
  //  1xg gxp pxr
  // -->  first is element-wise row multiplication
  //  results in a g x r matrix

  ArrayXd x1( ng );
  for (int i=0;i<ng;i++)
    x1[i] = sqrt( ( n * prior[i] ) * fac );

  
  // sqrt((n * prior)*fac) * scale(group.means, center = xbar, scale = FALSE)
  Eigen::MatrixXd group_means_centered( ng , p );
  for (int i=0;i<ng;i++)
    for (int j=0;j<p;j++)
      group_means_centered(i,j) = x1[i] * ( group_means(i,j) - xbar[j] ) ;
  
  Eigen::MatrixXd X2 = Eigen::MatrixXd::Zero( ng , rank );
  for (int i=0;i<ng;i++)
    for (int j=0;j<rank;j++)
      for (int k=0;k<p;k++)
	X2(i,j) += group_means_centered(i,k) * scaling2(k,j); 

  
  //
  // SVD of X2, g x r matrix
  //


  //okay = Statistics::svdcmp( X2 , W , V );     
  
  Eigen::BDCSVD<Eigen::MatrixXd> svd2( X2 , Eigen::ComputeThinU | Eigen::ComputeThinV );

  V = svd2.matrixV();
  W = svd2.singularValues();
 
  // if ( ! okay ) 
  //   {
  //     std::cout << "group means\n"
  // 		<< group_means.print() << "\n";
  //     std::cout << "X\n" << X2.print() << "\n";
  //     Helper::halt( "problem in lda_t, SVD 2" );
  //   }

  int rank2 = 0;
  for (int j=0;j<W.size();j++) 
    if ( W[j] > ( tol * W[0] ) ) ++rank2;
  
  if ( rank2 == 0 ) {
    model.valid = false;
    model.errmsg = "group means are numerically identical" ;
    return model;
  }
  
  //
  // ( p x r ) . ( r , r2 ) 
  // scaling <- scaling %*% X.s$v[, 1L:rank]
  //

  Eigen::MatrixXd scaling3 = Eigen::MatrixXd::Zero( p , rank2 );
  for (int i=0;i<p;i++)
    for (int j=0;j<rank2;j++)
      for (int k=0;k<rank;k++)
	scaling3(i,j) += scaling2(i,k) * V(k,j);
  
  model.valid = true;
  model.prior = prior;
  model.counts = counts;
  model.means = group_means;
  model.scaling = scaling3;
  model.n = n;
  model.svd.resize( rank2 );
  for (int i=0;i<rank2;i++) 
    model.svd[i] = W[i];

  model.labels.clear();
  
  std::map<std::string,int>::const_iterator ll = counts.begin();
  while ( ll != counts.end() )
    {
      model.labels.push_back( ll->first );
      ++ll;
    }
  
  return model;
}
  

lda_posteriors_t lda_t::predict( const lda_model_t & model , const Eigen::MatrixXd & X )
{

  const int p = X.cols();
  const int n = X.rows();

  //  std::cout << "p = " << p << " " << n << " " << model.means.cols() << "\n";
  
  if ( p != model.means.cols() )
    Helper::halt( "wrong number of columns in lda_t::predict()" );  

  // use prior from training/model
  const int ng = model.prior.size();
  
  // remove overall means to keep distances small
  ArrayXd means = Eigen::ArrayXd::Zero( p );

  for (int i=0;i<ng;i++)
    for (int j=0;j<p;j++)
      means[j] += model.prior[i] * model.means(i,j) ;
  
  Eigen::MatrixXd X1(n,p); 
  for (int i=0;i<n;i++)
    for (int j=0;j<p;j++)
      X1(i,j) = X(i,j) - means(j);

  // n.p X p.r --> n.r
  // x <- scale(x, center = means, scale = FALSE) %*% scaling
  Eigen::MatrixXd X2 = X1 * model.scaling;
    
  Eigen::MatrixXd M1(ng,p); 
  for (int i=0;i<ng;i++)
    for (int j=0;j<p;j++)
      M1(i,j) = model.means(i,j) - means(j);
  
  
  // dm <- scale(object$means, center = means, scale = FALSE) %*% scaling
  Eigen::MatrixXd DM = M1 * model.scaling;

  // can be a lower dimension, of first N components if needed
  // for now, fix to N components
  // dimen <- if(missing(dimen)) length(object$svd) else min(dimen, length(object$svd))
  //     dm <- dm[, 1L:dimen, drop = FALSE]
  // for now, we keep all dimensions, and so no changes to DM

  const int dimen = model.svd.size();

  //     dist <- matrix(0.5 * rowSums(dm^2) - log(prior), nrow(x),
  //                    length(prior), byrow = TRUE) - x[, 1L:dimen, drop=FALSE] %*% t(dm)

  //   ( NxNG ) - ( n x dimen ) . ( dimen x NG   )   
  // --> NxNG

  Eigen::MatrixXd dist( n , ng );
  // priors
  Data::Vector<double> rowSums( ng );
  
  for (int i=0;i<ng;i++)
    {
      rowSums[i] = 0;
      for (int k=0;k<dimen;k++) 
	rowSums[i] += DM(i,k) * DM(i,k);
    }
  
  Data::Vector<double> p0( ng );
  for (int i=0;i<ng;i++) 
    p0[i] = 0.5 * rowSums[i] - log( model.prior[i] );

  for (int i=0;i<n;i++)
    for (int j=0;j<ng;j++)
      dist(i,j) = p0[j];
  
  for (int i=0;i<n;i++)
    for (int j=0;j<ng;j++)
      for (int k=0;k<dimen;k++)
	dist(i,j) -= X2(i,k) * DM( j , k ) ; // io.e. t(DM)
  
  //     dist <- exp( -(dist - apply(dist, 1L, min, na.rm=TRUE)))
  
  // normalize to 0 each indiv.
  for (int i=0;i<n;i++)
    {
      double mn = dist(i,0);
      for (int j=1;j<ng;j++) if ( dist(i,j) < mn ) mn = dist(i,j) ;
      double s = 0;
      for (int j=0;j<ng;j++) { dist(i,j) = exp( - ( dist(i,j) - mn ) ); s += dist(i,j); }
      for (int j=0;j<ng;j++) dist(i,j) /= s;
    }
  
  std::vector<int> cli( n );
  std::vector<std::string> cl( n );
  for (int i=0;i<n;i++)
    {
      int mxi = 0;
      for (int k=1;k<ng;k++) if ( dist(i,k) > dist(i,mxi) ) mxi = k;
      cli[i] = mxi;
      cl[i] = model.labels[ mxi ];
    }
  

  //
  // Construct final return object
  //
  			      
  lda_posteriors_t pp;

  // posterior probs
  pp.pp = dist;
  
  // most likely class assignment (idx)
  pp.cli = cli;
  
  // as above, but using string label from model
  pp.cl = cl;
  
  return pp;

} 
