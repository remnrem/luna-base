#include <stdio.h>
#include <stdlib.h>
#include "libICA.h"

/*
 * Creates a matrix
 */
static double **mat_create(int rows, int cols)
{
	double **M; int i;
	
	M = (double**) malloc(rows * sizeof(double*));
	for (i=0; i<rows; i++)
		M[i] = (double*) malloc(cols * sizeof(double));
	
	return M;
}

/*
 * Reads matrix M from fp
 */
static double **mat_read(FILE *fp, int *rows, int *cols)
{
	int i, j; double **M;

	fscanf(fp, "%d %d", rows, cols);
	M = mat_create(*rows, *cols);
	
	for (i=0; i<*rows; i++) {
		for (j=0; j<*cols; j++)
			fscanf(fp, "%lf ", &(M[i][j]));	
	}
	
	return M;	
}

/*
 * Prints matrix M to stdout
 */
static void mat_print(double **M, int rows, int cols)
{
	int i, j;

	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++) {
			printf("%0.6f", M[i][j]);
			if (j < cols - 1)
				printf(" ");
		}
		printf("\n");
	}
}


/*
 * Main program
 */
int main(int argc, char *argv[])
{
	double **X, **K, **W, **A, **S;
	int rows, cols, compc;
	FILE *fp;

	// Input parameters check
	if (argc == 2) {
		if ((fp = fopen(argv[1], "r")) == NULL) {
			perror("Error opening input file");
			exit(-1);
		}
	} else {	
		printf("usage: %s data_file\n", argv[0]);
		exit(-1);
	}

	compc = 2;
	
	// Matrix creation
	X = mat_read(fp, &rows, &cols);
	W = mat_create(compc, compc);
	A = mat_create(compc, compc);
	K = mat_create(cols, compc);
	S = mat_create(rows, cols);	
	
	// ICA computation
	fastICA(X, rows, cols, compc, K, W, A, S);

	// Output
	printf("$K\n");	
	mat_print(K, cols, compc);

	printf("\n$W\n");
	mat_print(W, compc, compc);

	printf("\n$A\n");
	mat_print(A, compc, compc);
	
	printf("\n$S\n");
	mat_print(S, rows, compc);
	
	return 0;
}
