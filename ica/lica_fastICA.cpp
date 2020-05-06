/**
 * @file libICA.c
 * 
 * Main FastICA functions.
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "ica/ica.h"
#include "ica/lica_fastICA.h"
#include "ica/lica_matrix.h"

#include "helper/logger.h"

extern logger_t logger;

void print( Data::Matrix<double> & D , int r = 0 , int c = 0 )
{
  if ( r == 0 ) r = D.dim1();
  if ( c == 0 ) c = D.dim2();

  printf("\n");
  for (int i=0;i<r;i++)
    {
      for (int j=0;j<c;j++)
	printf( " %f" , D[i][j] );
      printf("\n");
    }
}

void printv( Data::Vector<double> & D , int r = 0  )
{
  if ( r == 0 ) r = D.size();

  printf("\n");
  for (int i=0;i<r;i++)
    printf( " %f" , D[i] );
  printf("\n");
}



/*
 * Functions for matrix elements transformations
 * used in mat_apply_fx().
 */
static scal fx_inv(scal x, scal par)
{
  return (1/x);
}

static scal fx_inv_sqrt(scal x, scal par)
{
  return (1/sqrt(x));	
}

static scal fx_div_c(scal x, scal par)
{
  return (x/par);
}

static scal fx_rand(scal x, scal par)
{
  return (scal)rand()/RAND_MAX; 
}

static scal fx_tanh(scal x, scal par)
{
  return tanh(LIBICA_ALPHA * x);
}

static scal fx_1sub_sqr(scal x, scal par)
{
  return (LIBICA_ALPHA * (1-x*x));
}


/**
 * ICA function. Computes the W matrix from the
 * preprocessed data.
 */

static Data::Matrix<double> ICA_compute( Data::Matrix<double> & X )
{

  const int rows = X.dim1();
  const int cols = X.dim2();
  
  //  mat TXp, GWX, W, Wd, W1, D, TU, TMP;
  //vect d, lim;
  //int i, it;
  
  printf( "%i rows x %i cols \n" , rows , cols );

  // matrix creation
//   TXp = mat_create(cols, rows);
//   GWX = mat_create(rows, cols);
//   W = mat_create(rows, rows);
//   Wd = mat_create(rows, rows);
//   D = mat_create(rows, rows);
//   TMP = mat_create(rows, rows);
//   TU = mat_create(rows, rows);
//   W1 = mat_create(rows, rows);
//   d = vect_create(rows);
  

  Data::Matrix<double> TXp(cols, rows);
  Data::Matrix<double> GWX(rows, cols);
  Data::Matrix<double> W(rows, rows);
  Data::Matrix<double> Wd(rows, rows);
  Data::Matrix<double> D(rows, rows);
  Data::Matrix<double> TMP(rows, rows);
  Data::Matrix<double> TU(rows, rows);
  Data::Matrix<double> W1(rows, rows);
  Data::Vector<double> d( rows );

  // W rand init
  mat_apply_fx(W, fx_rand, 0);

  // sW <- La.svd(W)
  Wd = W;
  ica_t::cpp_svdcmp(Wd, d, D);
  
  // W <- sW$u %*% diag(1/sW$d) %*% t(sW$u) %*% W
  mat_transpose(Wd, TU);
  vect_apply_fx(d, fx_inv, 0);
  mat_diag(d, D);
  mat_mult(Wd, D, TMP);
  mat_mult(TMP, TU, D);
  mat_mult(D, W, Wd); // W = Wd

  // W1 <- W 
  W1 = Wd; // mat_copy(Wd, W1);
  
  // lim <- rep(1000, maxit); it = 1
  Data::Vector<double> lim( MAX_ITERATIONS );
  for (int i=0; i<MAX_ITERATIONS; i++)
    lim[i] = 1000;
  int it = 0;
  
  
  // t(X)/p
  mat_transpose(X, TXp);
  mat_apply_fx(TXp, fx_div_c, cols);

  while (lim[it] > TOLERANCE && it < MAX_ITERATIONS) {

    if ( it+1 % 50 == 0 ) 
      logger << " iteration " << it << " (f = " << lim[it] << ")\n";
    else
      logger << ".";
    
    // wx <- W %*% X
    mat_mult(Wd, X, GWX);
    		
    // gwx <- tanh(alpha * wx)
    mat_apply_fx(GWX, fx_tanh, 0);

    // v1 <- gwx %*% t(X)/p
    mat_mult(GWX, TXp, TMP); // V1 = TMP
		
    // g.wx <- alpha * (1 - (gwx)^2)
    mat_apply_fx(GWX, fx_1sub_sqr, 0);

        
    // v2 <- diag(apply(g.wx, 1, FUN = mean)) %*% W
    mat_mean_rows(GWX,  d);
    mat_diag(d, D);
    mat_mult(D, Wd, TU); // V2 = TU
    
    // W1 <- v1 - v2
    mat_sub(TMP, TU, W1);

    // sW1 <- La.svd(W1)
    W = W1; //mat_copy(W1, rows, rows, W);
    ica_t::cpp_svdcmp(W, d, D); // sets W as U

    // W1 <- sW1$u %*% diag(1/sW1$d) %*% t(sW1$u) %*% W1
    mat_transpose(W, TU);
    vect_apply_fx(d, fx_inv, 0);
    mat_diag(d, D);
    //    print(D,2,2);
    mat_mult(W, D, TMP);
    mat_mult(TMP, TU, D);
    mat_mult(D, W1, W); // W1 = W

    // lim[it + 1] <- max(Mod(Mod(diag(W1 %*% t(W))) - 1))
    mat_transpose(Wd, TU);
    mat_mult(W, TU, TMP);

    lim[it+1] = fabs(mat_max_abs_diag(TMP) - 1);

    // W <- W1
    Wd = W; // mat_copy(W, rows, rows, Wd);
    
    it++;

  }

  return Wd;
}


/**
 * Main FastICA function. Centers and whitens the input
 * matrix, calls the ICA computation function ICA_compute()
 * and computes the output matrixes.
 */
void fastICA( Data::Matrix<double> & X , 
	      int compc, 
	      Data::Matrix<double> & K , 
	      Data::Matrix<double> & W , 
	      Data::Matrix<double> & A , 
	      Data::Matrix<double> & S )
{

  const int rows = X.dim1();
  const int cols = X.dim2();

  Data::Matrix<double> XT( cols, rows );
  Data::Matrix<double> V( cols , cols );
  Data::Matrix<double> TU( cols , cols ); 
  Data::Matrix<double> D( cols , cols );
  Data::Matrix<double> X1( compc , rows );

  Data::Vector<double> scale( cols );
  Data::Vector<double> d( cols );

  // matrix creation
//   XT = mat_create(cols, rows);
//   X1 = mat_create(compc, rows);
//   V = mat_create(cols, cols);
//   D = mat_create(cols, cols);
//   TU = mat_create(cols, cols);

//   scale = vect_create(cols);
//   d = vect_create(cols);


  printf( "  pre-processing...\n" );	

  /*
   * CENTERING
   */

  mat_center(X, scale);
  
  /*
   * WHITENING
   */
  
  // X <- t(X); V <- X %*% t(X)/rows 
  mat_transpose(X, XT);
  mat_apply_fx(X, fx_div_c, rows);
  mat_mult(XT, X, V);


  // La.svd(V)
  ica_t::cpp_svdcmp(V, d, D);  // V = s$u, d = s$d, D = s$v

  // D <- diag(c(1/sqrt(d))
  vect_apply_fx(d, fx_inv_sqrt, 0);	
  mat_diag(d, D);

  // K <- D %*% t(U)
  mat_transpose(V, TU);

  mat_mult(D, TU, V); // K = V 

  // take only first compc components
  Data::Matrix<double> V1( compc , cols );
  for (int r=0;r<compc;r++)
    for (int c=0;c<cols;c++)
      V1[r][c] = V[r][c];
   
  // X1 <- K %*% X
  mat_mult(V1, XT, X1);


  /*
   * FAST ICA
   */
  
  logger << "    starting ICA\n";
  
  Data::Matrix<double> _A = ICA_compute(X1);
  

  /*
   * OUTPUT
   */

  // X <- t(x)
  mat_transpose(XT, X);
  mat_decenter(X, scale);

  // K
  mat_transpose(V1, K);

  // w <- a %*% K; S <- w %*% X
  D.resize( compc , cols );
  mat_mult(_A, V1, D);	
  mat_mult(D, XT, X1);

  // S
  mat_transpose(X1, S);

  // A <- t(w) %*% solve(w * t(w))
  D.resize(compc,compc);
  TU.resize(compc,compc);
  mat_transpose(D, TU);

  Data::Matrix<double> V2( compc , compc );
  mat_mult(D, TU, V2);
  mat_inverse(V2, D);
  mat_mult(TU, D, V2);

  // A
  mat_transpose(V2, A);

  // W
  mat_transpose(_A, W);

}
