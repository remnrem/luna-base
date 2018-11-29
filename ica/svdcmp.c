/** 
 * @file svdcmp.c
 * 
 * Singular value decomposition.
 */

#include <stdio.h>
#include <math.h>
#include "matrix.h"

#define true  1
#define false 0


/**
 * Computes the singular value decomposition of A.
 * 
 * Give a matrix A, dimensions M by N, this routine computes its
 * singular value decomposition, A = U * W * transpose V. The matrix U
 * replaces A on output. The diagonal matrix of singular values, W, is
 * output as a vector W. The matrix V (not the transpose of V) is output
 * as V. M must be greater or equal to N. If it is smaller then A should
 * be filled up to square with zero rows.
 */
void svdcmp(mat A, int M, int N, vect W, mat V)
{
    /* Householder reduction to bidiagonal form. */
    int NM;
    scal C;
    scal F;
    scal G = 0.0;
    scal H;
    scal S;
    scal X;
    scal Y;
    scal Z;
    scal Scale = 0.0;
    scal ANorm = 0.0;
    scal tmp;
    int flag;
    int i;
    int its;
    int j;
    int jj;
    int k;
    int l;
	vect rv1;


	rv1 = vect_create(N);

    if( M < N ) {
        fprintf( stderr, "You must augment A with extra zero rows.\n" );
        return;
    }

    for( i = 0; i < N; ++i ) {
        l = i + 1;
        rv1[i] = Scale * G;
        G = 0.0;
        S = 0.0;
        Scale = 0.0;
        if( i < M ) {
            for( k = i; k < M; ++k ) {
                Scale = Scale + fabs( A[k][i] );
            }
            if( Scale != 0.0 ) {
                for( k = i; k < M; ++k ) {
                    A[k][i] = A[k][i] / Scale;
                    S = S + A[k][i] * A[k][i];
                }
                F = A[i][i];
                G = sqrt(S);
                if( F > 0.0 ) {
                    G = -G;
                }
                H = F * G - S;
                A[i][i] = F - G;
                if( i != (N-1) ) {
                    for( j = l; j < N; ++j ) {
                        S = 0.0;
                        for( k = i; k < M; ++k ) {
                            S = S + A[k][i] * A[k][j];
                        }
                        F = S / H;
                        for( k = i; k < M; ++k ) {
                            A[k][j] = A[k][j] + F * A[k][i];
                        }
                    }
                }
                for( k = i; k < M; ++k ) {
                    A[k][i] = Scale * A[k][i];
                }
            }
        }

        W[i] = Scale * G;
        G = 0.0;
        S = 0.0;
        Scale = 0.0;
        if( (i < M) && (i != (N-1)) ) {
            for( k = l; k < N; ++k ) {
                Scale = Scale + fabs( A[i][k] );
            }
            if( Scale != 0.0 ) {
                for( k = l; k < N; ++k ) {
                    A[i][k] = A[i][k] / Scale;
                    S = S + A[i][k] * A[i][k];
                }
                F = A[i][l];
                G = sqrt(S);
                if( F > 0.0 ) {
                    G = -G;
                }
                H = F * G - S;
                A[i][l] = F - G;
                for( k = l; k < N; ++k ) {
                    rv1[k] = A[i][k] / H;
                }
                if( i != (M-1) ) {
                    for( j = l; j < M; ++j ) {
                        S = 0.0;
                        for( k = l; k < N; ++k ) {
                            S = S + A[j][k] * A[i][k];
                        }
                        for( k = l; k < N; ++k ) {
                            A[j][k] = A[j][k] + S * rv1[k];
                        }
                    }
                }
                for( k = l; k < N; ++k ) {
                    A[i][k] = Scale * A[i][k];
                }
            }
        }
        tmp = fabs( W[i] ) + fabs( rv1[i] );
        if( tmp > ANorm )
            ANorm = tmp;
    }

    /* Accumulation of right-hand transformations. */
    for( i = N-1; i >= 0; --i ) {
        if( i < (N-1) ) {
            if( G != 0.0 ) {
                for( j = l; j < N; ++j ) {
                    V[j][i] = (A[i][j] / A[i][l]) / G;
                }
                for( j = l; j < N; ++j ) {
                    S = 0.0;
                    for( k = l; k < N; ++k ) {
                        S = S + A[i][k] * V[k][j];
                    }
                    for( k = l; k < N; ++k ) {
                        V[k][j] = V[k][j] + S * V[k][i];
                    }
                }
            }
            for( j = l; j < N; ++j ) {
                V[i][j] = 0.0;
                V[j][i] = 0.0;
            }
        }
        V[i][i] = 1.0;
        G = rv1[i];
        l = i;
    }

    /* Accumulation of left-hand transformations. */
    for( i = N-1; i >= 0; --i ) {
        l = i + 1;
        G = W[i];
        if( i < (N-1) ) {
            for( j = l; j < N; ++j ) {
                A[i][j] = 0.0;
            }
        }
        if( G != 0.0 ) {
            G = 1.0 / G;
            if( i != (N-1) ) {
                for( j = l; j < N; ++j ) {
                    S = 0.0;
                    for( k = l; k < M; ++k ) {
                        S = S + A[k][i] * A[k][j];
                    }
                    F = (S / A[i][i]) * G;
                    for( k = i; k < M; ++k ) {
                        A[k][j] = A[k][j] + F * A[k][i];
                    }
                }
            }
            for( j = i; j < M; ++j ) {
                A[j][i] = A[j][i] * G;
            }
        } else {
            for( j = i; j < M; ++j ) {
                A[j][i] = 0.0;
            }
        }
        A[i][i] = A[i][i] + 1.0;
    }

    /* Diagonalization of the bidiagonal form.
       Loop over singular values. */
    for( k = (N-1); k >= 0; --k ) {
        /* Loop over allowed iterations. */
        for( its = 1; its <= 30; ++its ) {
            /* Test for splitting.
               Note that rv1[0] is always zero. */
            flag = true;
            for( l = k; l >= 0; --l ) {
                NM = l - 1;
                if( (fabs(rv1[l]) + ANorm) == ANorm ) {
                    flag = false;
                    break;
                } else if( (fabs(W[NM]) + ANorm) == ANorm ) {
                    break;
                }
            }

            /* Cancellation of rv1[l], if l > 0; */
            if( flag ) {
                C = 0.0;
                S = 1.0;
                for( i = l; i <= k; ++i ) {
                    F = S * rv1[i];
                    if( (fabs(F) + ANorm) != ANorm ) {
                        G = W[i];
                        H = sqrt( F * F + G * G );
                        W[i] = H;
                        H = 1.0 / H;
                        C = ( G * H );
                        S = -( F * H );
                        for( j = 0; j < M; ++j ) {
                            Y = A[j][NM];
                            Z = A[j][i];
                            A[j][NM] = (Y * C) + (Z * S);
                            A[j][i] = -(Y * S) + (Z * C);
                        }
                    }
                }
            }
            Z = W[k];
            /* Convergence. */
            if( l == k ) {
                /* Singular value is made nonnegative. */
                if( Z < 0.0 ) {
                    W[k] = -Z;
                    for( j = 0; j < N; ++j ) {
                        V[j][k] = -V[j][k];
                    }
                }
                break;
            }

            if( its >= 30 ) {
                fprintf( stderr, "No convergence in 30 iterations.\n" );
                return;
            }

            X = W[l];
            NM = k - 1;
            Y = W[NM];
            G = rv1[NM];
            H = rv1[k];
            F = ((Y-Z)*(Y+Z) + (G-H)*(G+H)) / (2.0*H*Y);
            G = sqrt( F * F + 1.0 );
            tmp = G;
            if( F < 0.0 )
                tmp = -tmp;
            F = ((X-Z)*(X+Z) + H*((Y/(F+tmp))-H)) / X;

            /* Next QR transformation. */
            C = 1.0;
            S = 1.0;
            for( j = l; j <= NM; ++j ) {
                i = j + 1;
                G = rv1[i];
                Y = W[i];
                H = S * G;
                G = C * G;
                Z = sqrt( F * F + H * H );
                rv1[j] = Z;
                C = F / Z;
                S = H / Z;
                F = (X * C) + (G * S);
                G = -(X * S) + (G * C);
                H = Y * S;
                Y = Y * C;
                for( jj = 0; jj < N; ++jj ) {
                    X = V[jj][j];
                    Z = V[jj][i];
                    V[jj][j] = (X * C) + (Z * S);
                    V[jj][i] = -(X * S) + (Z * C);
                }
                Z = sqrt( F * F + H * H );
                W[j] = Z;

                /* Rotation can be arbitrary if Z = 0. */
                if( Z != 0.0 ) {
                    Z = 1.0 / Z;
                    C = F * Z;
                    S = H * Z;
                }
                F = (C * G) + (S * Y);
                X = -(S * G) + (C * Y);
                for( jj = 0; jj < M; ++jj ) {
                    Y = A[jj][j];
                    Z = A[jj][i];
                    A[jj][j] = (Y * C) + (Z * S);
                    A[jj][i] = -(Y * S) + (Z * C);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = F;
            W[k] = X;
        }
    }
	vect_delete(rv1);
	
    return;
}
