
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


#include "ica/ica.h"

#include "stats/eigen_ops.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "eval.h"
#include "db/db.h"

extern writer_t writer;

extern logger_t logger;


bool eigen_ica_t::proc( Eigen::MatrixXd & X , int nc )
{

  int rows = X.rows();
  int cols = X.cols();
  
  if ( rows < 2 || cols < 2 ) 
    return false;
   
  W.resize( nc, nc );
  A.resize( nc, nc );
  K.resize( cols , nc );
  S.resize( rows, cols );
  
  fastICA( X, nc, W, A, K, S);
  
  return true;
}





void eigen_ica_t::fastICA( Eigen::MatrixXd & X , 
			   int nc , 
			   Eigen::MatrixXd & W , 
			   Eigen::MatrixXd & A , 
			   Eigen::MatrixXd & K , 
			   Eigen::MatrixXd & S ) 
{

  // W  nc   x  nc
  // A  nc   x  nc
  // K  cols x  nc 
  // S  rows x  cols
  
  // fun = logcosh
  // alpha = 1 

  // row.norm = F

  // inputs = X & n.comp

  const int n = X.rows();
  const int p = X.cols();

  const int minnp = n < p ? n : p;

  if ( nc > minnp )
    {
      logger << "  ** warning: nc is too large, resetting to " << minnp << "\n";
      nc = minnp;
    }
  
  

  //
  // Centering
  //
  
  eigen_ops::scale( X , row_norm );
  
  //
  // transpose
  //
  
  X.transposeInPlace();


  //
  // Whitening
  //

  // X %*% t(X)/n
  Eigen::MatrixXd V = X * ( X.array() / n ).matrix().transpose();

  // s <- La.svd(V)
  Eigen::BDCSVD<Eigen::MatrixXd> s( V , Eigen::ComputeThinU | Eigen::ComputeThinV );

  // U = sW.matrixU();
  // V = sW.matrixV();
  // W = sW.singularValues();

    // D <- diag(c(1/sqrt(s$d)))
  Eigen::MatrixXd D = s.singularValues().array().sqrt().inverse().matrix().asDiagonal();
    
  // K <- D %*% t(s$u)
  K = D * s.matrixU().transpose();
  
  // K <- matrix( K[1:n.comp, ], n.comp, p)	      
  K = K.block( 0 , 0 , nc , p );
    
  // X1 <- K %*% X
  Eigen::MatrixXd X1 = K * X;
  
  //
  // parallel method
  //

  // a <-  ica.R.par(X1, n.comp, tol = tol, fun = fun, alpha = alpha, 
  // 		  maxit = maxit, verbose = verbose, w.init = w.init)
  
  Eigen::MatrixXd a = ica_parallel( X1 , nc );

  //
  // Get results to return
  //
  
  // w <- a %*% K
  Eigen::MatrixXd w = a * K;

  // S <- w %*% X
  S = w * X;
  
  // A <- t(w) %*% solve(w %*% t(w))
  A = w.transpose() * ( w * w.transpose() ).inverse();
  
  // return(list(X = t(X), K = t(K), W = t(a), A = t(A), S = t(S)))

  K.transposeInPlace();
  W = a.transpose();
  A.transposeInPlace();
  S.transposeInPlace();

  std::cout << "K dims = " << K.rows() << " " << K.cols() << "\n";
  std::cout << "W dims = " << W.rows() << " " << W.cols() << "\n";
  std::cout << "A dims = " << A.rows() << " " << A.cols() << "\n";
  std::cout << "S dims = " << S.rows() << " " << S.cols() << "\n";
  
  std::cerr << " all done\n";
}




//
// ICA PAR
//


Eigen::MatrixXd eigen_ica_t::ica_parallel( const Eigen::MatrixXd & X ,
 					   const int nc )
{
  
	 				   
  const int p = X.cols();
  
  //
  // initialize W with random normal values
  //
  
  Eigen::MatrixXd W( nc , nc );

  eigen_ops::random_normal( W );

  //
  // SVD
  //

  //  sW <- La.svd(W)

  Eigen::BDCSVD<Eigen::MatrixXd> sW( W , Eigen::ComputeThinU | Eigen::ComputeThinV );

  // W <- sW$u %*% Diag(1/sW$d) %*% t(sW$u) %*% W
  W = sW.matrixU() * ( sW.singularValues().array().inverse() ).matrix().asDiagonal() * sW.matrixU().transpose() * W;
  
  Eigen::MatrixXd W1 = W;

  //lim <- rep(1000, maxit)
  std::vector<double> lim( maxit , 1000 );
  
  // iteration counter
  int it = 0;

  logger << "  starting iteraatons (symmetric FastICA using logcosh approx. to neg-entropy function)";
  
  while ( lim[it] > tol && it < (maxit-1) )
    {
      //  wx <- W %*% X
      // gwx <- tanh(alpha * wx)
      // alpha = 1 , so ignore

      Eigen::MatrixXd gwx = ( W * X ).array().tanh().matrix();

      // v1 <- gwx %*% t(X)/p
      Eigen::MatrixXd v1 = gwx * ( X.array() / p ).matrix().transpose(); 
      
      // g.wx <- alpha * (1 - (gwx)^2)
      // nb alpha == 1
      gwx = 1 - gwx.array().square(); 

      //v2 <- Diag(apply(g.wx, 1, FUN = mean)) %*% W
      Eigen::MatrixXd v2 = gwx.array().rowwise().mean().matrix().asDiagonal() * W;
            
      //W1 <- v1 - v2
      W1 = v1 - v2;
      
      // sW1 <- La.svd(W1)      
      Eigen::BDCSVD<Eigen::MatrixXd> sW1( W1 , Eigen::ComputeThinU | Eigen::ComputeThinV );
      
      // W1 <- sW1$u %*% Diag(1/sW1$d) %*% t(sW1$u) %*% W1
      W1 = sW1.matrixU() * ( sW1.singularValues().array().inverse() ).matrix().asDiagonal() * sW1.matrixU().transpose() * W1;
      
      // lim[it + 1] <- max( Mod(   Mod(  diag(W1 %*% t(W) )  )  - 1 ) )
      lim[ it + 1 ] =  ( ( W1 * W.transpose() ).diagonal().array().abs() - 1 ).abs().maxCoeff();

      // W <- W1
      W = W1;

      if ( it % 50 == 0 ) logger << "\n ";
      if ( it % 10 == 0 ) logger << " ";
      logger << ".";
      
      // next iteration
      ++it;
    }

    //
    // Return W
    //

    return W;

}


