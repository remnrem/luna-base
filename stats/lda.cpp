
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


lda_model_t lda_t::fit()
{
  lda_model_t model;

  const int n = X.dim1();
  const int p = X.dim2();

  std::map<std::string,int> counts;
  for (int i=0;i<y.size();i++) counts[ y[i] ]++;

  // assign index to groups
  std::map<std::string,int> gidx;
  std::map<std::string,int>::const_iterator cc = counts.begin();
  while ( cc != counts.end() )
    {
      //std::cout << " cc -> " << cc->first << " " << cc->second << "\n";
      const int sz = gidx.size();
      gidx[ cc->first ] = sz;
      ++cc;
    }

  std::vector<int> yi( n );
  for (int i=0;i<n;i++) { 
    yi[i] = gidx[y[i]];    
  }
    
  // number of classes
  const int ng = counts.size();
  
  // priors
  std::vector<double> prior;
  cc = counts.begin();
  while ( cc != counts.end() )
    {
      prior.push_back( cc->second / (double)n );
      ++cc;
    }


  // group means
  Data::Matrix<double> group_means( ng , p );
  for ( int i = 0 ; i < n ; i++ )
    for ( int j = 0 ; j < p ; j++ )
      group_means( yi[i], j ) += X(i,j);

  for ( int i = 0 ; i < ng ; i++ )
    for ( int j = 0 ; j < p ; j++ )
      group_means( i, j ) /= n * prior[i];

  //  std::cout << "GM\n" << group_means.print() << "\n";
  
  // adjust X by group mean; get variance of each measure (i.e. looking for within-group variability)
  std::vector<double> f1(p);
  for (int j=0;j<p;j++)
    {
      std::vector<double> xa(n);
      for (int i=0;i<n;i++) xa[i] = X(i,j) - group_means( yi[i] , j ) ;
      f1[j] = MiscMath::sdev( xa );      
      if ( f1[j] < tol ) {
	model.valid = false;
	model.errmsg = "variable " + Helper::int2str(j) + " is constant within group ";
	return model;
      }      
    }

  
  // scaling matrix (diagonal)... can keep as vector
  Data::Matrix<double> scaling( p , p , 0 );
  for (int j=0;j<p;j++) scaling(j,j) = 1.0 / f1[j] ;
  
  double fac = 1.0 / double( n - ng );
  double sqrt_fac = sqrt( fac );
  
  // scaled X
  Data::Matrix<double> X1( n , p );
  for (int i=0;i<n;i++)
    for (int j=0;j<p;j++)      
      X1(i,j) = sqrt_fac * ( X(i,j) - group_means( yi[i] , j ) ) * scaling(j,j); 
  

  Data::Vector<double> W(p);
  Data::Matrix<double> V(p,p);
  
  bool okay = Statistics::svdcmp( X1 , W , V );

  // SVD (no left singular vectors needed)
  // X.s <- svd(X1, nu = 0L)

  // order so that we have descending singular values

  std::vector<int> o;
  std::vector<bool> in( W.size() , false ) ;
  for (int i=0;i<W.size();i++) 
    {
      int mxi = 0;
      for (int j=0;j<W.size();j++) { if ( ! in[j] ) { mxi = j ; break; } }
      for (int j=0;j<W.size();j++) if ( (!in[j]) && W[j] >= W[mxi] ) mxi = j ;      
      o.push_back( mxi );
      in[ mxi ] = true;
    }

  Data::Vector<double> W2 = W;
  Data::Matrix<double> V2 = V;
  for (int i=0;i<W.size();i++)
    W[i] = W2[ o[i] ];
  
  // swap cols of V
  for (int i=0;i<p;i++)
    for (int j=0;j<p;j++)
      V(i,j) = V2( i , o[j] );
  
  // rank
  int rank = 0;
  for (int j=0;j<p;j++) if ( W[j] > tol ) ++rank;

  if ( rank == 0 ) {
    model.valid = false;
    model.errmsg = "problem with collinearity/constant values in input data" ;
    return model;
  }

  if ( rank < p ) logger << " warning... rank < p\n"; 

  // get new scaling matrix p x rank , possibly of rank < p
  // scaling <- scaling %*% X.s$v[, 1L:rank] %*% diag(1/X.s$d[1L:rank],,rank)


    
  Data::Matrix<double> scaling2( p , rank );
  for (int i=0;i<p;i++)
    for (int j=0;j<rank;j++)
      for (int k=0;k<p;k++)
	scaling2(i,j) += scaling(i,k) * V(k,j);

  for (int i=0;i<p;i++)
    for (int j=0;j<rank;j++)
      scaling2(i,j) /= W[j];
	  
   // mean values of each measure
  Data::Vector<double> xbar( p );
  for (int i=0;i<ng;i++)
    for (int j=0;j<p;j++)
      xbar[j] += prior[i] * group_means(i,j);
  
  
  fac = 1.0/(double)(ng - 1);
  
  //X <- sqrt((n * prior)*fac) * scale(group.means, center = xbar, scale = FALSE) %*% scaling
  //  1xg gxp pxr
  // -->  first is element-wise row multiplication
  //  results in a g x r matrix

  Data::Vector<double> x1( ng );
  for (int i=0;i<ng;i++)
    x1[i] = sqrt( ( n * prior[i] ) * fac );

  // sqrt((n * prior)*fac) * scale(group.means, center = xbar, scale = FALSE)
  Data::Matrix<double> group_means_centered( ng , p );
  for (int i=0;i<ng;i++)
    for (int j=0;j<p;j++)
      group_means_centered(i,j) = x1[i] * ( group_means(i,j) - xbar[j] ) ;
  
  Data::Matrix<double> X2( ng , rank );
  for (int i=0;i<ng;i++)
    for (int j=0;j<rank;j++)
      for (int k=0;k<p;k++)
	X2(i,j) += group_means_centered(i,k) * scaling2(k,j); 
  

  // SVD of X2, g x r matrix
  
  W.clear();
  V.clear();
  
  W.resize( rank );
  V.resize( rank , rank );
        
  okay = Statistics::svdcmp( X2 , W , V );
  
  int rank2 = 0;
  for (int j=0;j<rank;j++) if ( W[j] > tol ) ++rank2;
  
  // re-order elements once more
  o.clear();
  in.clear();
  in.resize( W.size() , false ) ;
  for (int i=0;i<W.size();i++) 
    {
      int mxi = -1;
      for (int j=0;j<W.size();j++) { if ( ! in[j] ) { mxi = j ; break; } } 
      for (int j=0;j<W.size();j++) if ( (!in[j]) && W[j] >= W[mxi] ) mxi = j ;      
      o.push_back( mxi );
      in[ mxi ] = true;
    }

  W2 = W;
  V2 = V;

  for (int i=0;i<W.size();i++)
    W[i] = W2[ o[i] ]; 
  
  // for V
  for (int i=0;i<rank;i++)
    for (int j=0;j<rank;j++)
      V(i,j) = V2( i , o[j] );
    
  if ( rank2 == 0 ) {
    model.valid = false;
    model.errmsg = "group means are numerically identical" ;
    return model;
  }

  
  // ( p x r ) . ( r , r2 ) 
  // scaling <- scaling %*% X.s$v[, 1L:rank]

  Data::Matrix<double> scaling3( p , rank2 );
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
  for (int i=0;i<rank2;i++) model.svd[i] = W[i];

  model.labels.clear();
  
  std::map<std::string,int>::const_iterator ll = counts.begin();
  while ( ll != counts.end() )
    {
      model.labels.push_back( ll->first );
      ++ll;
    }

  return model;
}
  

lda_posteriors_t lda_t::predict( const lda_model_t & model , const Data::Matrix<double> & X )
{

  const int p = X.dim2();
  const int n = X.dim1();

  //  std::cout << "p = " << p << " " << n << " " << model.means.dim2() << "\n";
  
  if ( p != model.means.dim2() )
    Helper::halt( "wrong number of columns in lda_t::predict()" );  

  // use prior from training/model
  const int ng = model.prior.size();
  
  // remove overall means to keep distances small
  Data::Vector<double> means( p );
  for (int i=0;i<ng;i++)
    for (int j=0;j<p;j++)
      means[j] += model.prior[i] * model.means(i,j) ;

  Data::Matrix<double> X1(n,p); 
  for (int i=0;i<n;i++)
    for (int j=0;j<p;j++)
      X1(i,j) = X(i,j) - means(j);

  // n.p X p.r --> n.r
  // x <- scale(x, center = means, scale = FALSE) %*% scaling
  Data::Matrix<double> X2 = X1 * model.scaling;
    
  Data::Matrix<double> M1(ng,p); 
  for (int i=0;i<ng;i++)
    for (int j=0;j<p;j++)
      M1(i,j) = model.means(i,j) - means(j);
  
  
  // dm <- scale(object$means, center = means, scale = FALSE) %*% scaling
  Data::Matrix<double> DM = M1 * model.scaling;

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

  Data::Matrix<double> dist( n , ng );
  // priors
  Data::Vector<double> rowSums( ng );
  for (int i=0;i<ng;i++)
    for (int k=0;k<dimen;k++) rowSums[i] += DM(i,k) * DM(i,k);

  Data::Vector<double> p0( ng );
  for (int i=0;i<ng;i++) p0[i] = 0.5 * rowSums[i] - log( model.prior[i] );

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
  
  			      
  lda_posteriors_t pp;

  // posterior probs
  pp.pp = dist;

  // most likely class assignment (idx)
  pp.cli = cli;
  
  // as above, but using string label from model
  pp.cl = cl;
  
  return pp;

}
  

