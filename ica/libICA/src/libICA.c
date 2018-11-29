/**
 * @file libICA.c
 * 
 * Main FastICA functions.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libICA.h"
#include "matrix.h"
#include "svdcmp.h"


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
	return tanh(ALPHA * x);
}

static scal fx_1sub_sqr(scal x, scal par)
{
	return (ALPHA * (1-x*x));
}


/**
 * ICA function. Computes the W matrix from the
 * preprocessed data.
 */
static mat ICA_compute(mat X, int rows, int cols)
{
	mat TXp, GWX, W, Wd, W1, D, TU, TMP;
	vect d, lim;
	int i, it;

	// matrix creation
	TXp = mat_create(cols, rows);
	GWX = mat_create(rows, cols);
	W = mat_create(rows, rows);
	Wd = mat_create(rows, rows);
	D = mat_create(rows, rows);
	TMP = mat_create(rows, rows);
	TU = mat_create(rows, rows);
	W1 = mat_create(rows, rows);
	d = vect_create(rows);

	// W rand init
	mat_apply_fx(W, rows, rows, fx_rand, 0);

	// sW <- La.svd(W)
	mat_copy(W, rows, rows, Wd);
	svdcmp(Wd, rows, rows, d, D);

	// W <- sW$u %*% diag(1/sW$d) %*% t(sW$u) %*% W
	mat_transpose(Wd, rows, rows, TU);
	vect_apply_fx(d, rows, fx_inv, 0);
	mat_diag(d, rows, D);
	mat_mult(Wd, rows, rows, D, rows, rows, TMP);
	mat_mult(TMP, rows, rows, TU, rows, rows, D);
	mat_mult(D, rows, rows, W, rows, rows, Wd); // W = Wd

	// W1 <- W 
	mat_copy(Wd, rows, rows, W1);

	// lim <- rep(1000, maxit); it = 1
	lim = vect_create(MAX_ITERATIONS);
	for (i=0; i<MAX_ITERATIONS; i++)
		lim[i] = 1000;
	it = 0;


	// t(X)/p
	mat_transpose(X, rows, cols, TXp);
	mat_apply_fx(TXp, cols, rows, fx_div_c, cols);

	while (lim[it] > TOLERANCE && it < MAX_ITERATIONS) {
		// wx <- W %*% X
		mat_mult(Wd, rows, rows, X, rows, cols, GWX);

		// gwx <- tanh(alpha * wx)
		mat_apply_fx(GWX, rows, cols, fx_tanh, 0);
		
		// v1 <- gwx %*% t(X)/p
		mat_mult(GWX, rows, cols, TXp, cols, rows, TMP); // V1 = TMP
		
		// g.wx <- alpha * (1 - (gwx)^2)
		mat_apply_fx(GWX, rows, cols, fx_1sub_sqr, 0);

		// v2 <- diag(apply(g.wx, 1, FUN = mean)) %*% W
		mat_mean_rows(GWX, rows, cols, d);
		mat_diag(d, rows, D);
		mat_mult(D, rows, rows, Wd, rows, rows, TU); // V2 = TU

		// W1 <- v1 - v2
		mat_sub(TMP, TU, rows, rows, W1);
		
		// sW1 <- La.svd(W1)
		mat_copy(W1, rows, rows, W);
		svdcmp(W, rows, rows, d, D);

		// W1 <- sW1$u %*% diag(1/sW1$d) %*% t(sW1$u) %*% W1
		mat_transpose(W, rows, rows, TU);
		vect_apply_fx(d, rows, fx_inv, 0);
		mat_diag(d, rows, D);
		mat_mult(W, rows, rows, D, rows, rows, TMP);
		mat_mult(TMP, rows, rows, TU, rows, rows, D);
		mat_mult(D, rows, rows, W1, rows, rows, W); // W1 = W
		
		// lim[it + 1] <- max(Mod(Mod(diag(W1 %*% t(W))) - 1))
		mat_transpose(Wd, rows, rows, TU);
		mat_mult(W, rows, rows, TU, rows, rows, TMP);
		lim[it+1] = fabs(mat_max_diag(TMP, rows, rows) - 1);

		// W <- W1
		mat_copy(W, rows, rows, Wd);

		it++;
	}

	// clean up
	mat_delete(TXp, cols, rows);
	mat_delete(GWX, rows, cols);
	mat_delete(W, rows, rows);
	mat_delete(D, rows, rows);
	mat_delete(TMP, rows, rows);
	mat_delete(TU, rows, rows);
	mat_delete(W1, rows, rows);
	vect_delete(d);	

	return Wd;
}


/**
 * Main FastICA function. Centers and whitens the input
 * matrix, calls the ICA computation function ICA_compute()
 * and computes the output matrixes.
 */
void fastICA(mat X, int rows, int cols, int compc, mat K, mat W, mat A, mat S)
{
	mat XT, V, TU, D, X1, _A;
	vect scale, d;

	// matrix creation
	XT = mat_create(cols, rows);
	X1 = mat_create(compc, rows);
	V = mat_create(cols, cols);
	D = mat_create(cols, cols);
	TU = mat_create(cols, cols);
	scale = vect_create(cols);
	d = vect_create(cols);

	/*
	 * CENTERING
	 */
	mat_center(X, rows, cols, scale);


	/*
	 * WHITENING
	 */

	// X <- t(X); V <- X %*% t(X)/rows 
	mat_transpose(X, rows, cols, XT);
	mat_apply_fx(X, rows, cols, fx_div_c, rows);
	mat_mult(XT, cols, rows, X, rows, cols, V);
	
	// La.svd(V)
	svdcmp(V, cols, cols, d, D);  // V = s$u, d = s$d, D = s$v

	// D <- diag(c(1/sqrt(d))
	vect_apply_fx(d, cols, fx_inv_sqrt, 0);	
	mat_diag(d, cols, D);

	// K <- D %*% t(U)
	mat_transpose(V, cols, cols, TU);
	mat_mult(D, cols, cols, TU, cols, cols, V); // K = V 

	// X1 <- K %*% X
	mat_mult(V, compc, cols, XT, cols, rows, X1);

	/*
	 * FAST ICA
	 */
	_A = ICA_compute(X1, compc, rows);

	
	/*
	 * OUTPUT
	 */

	// X <- t(x)
	mat_transpose(XT, cols, rows, X);
	mat_decenter(X, rows, cols, scale);

	// K
	mat_transpose(V, compc, cols, K);

	// w <- a %*% K; S <- w %*% X
	mat_mult(_A, compc, compc, V, compc, cols, D);	
	mat_mult(D, compc, cols, XT, cols, rows, X1);
	// S
	mat_transpose(X1, compc, rows, S);

	// A <- t(w) %*% solve(w * t(w))
	mat_transpose(D, compc, compc, TU);
	mat_mult(D, compc, compc, TU, compc, compc, V);
	mat_inverse(V, compc, D);
	mat_mult(TU, compc, compc, D, compc, compc, V);
	// A
	mat_transpose(V, compc, compc, A);

	// W
	mat_transpose(_A, compc, compc, W);

	// cleanup
	mat_delete(XT, cols, rows);
	mat_delete(X1, compc, rows);
	mat_delete(V, cols, cols);
	mat_delete(D, cols, cols);
	mat_delete(TU,cols, cols);
	vect_delete(scale);
	vect_delete(d);
}
