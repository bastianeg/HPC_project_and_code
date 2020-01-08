#include <stdlib.h>
#include <stdio.h>
#include "matmult.h"

// allocate a double-prec m x n matrix
double ** malloc_2d(int m, int n) {
    if (m <= 0 || n <= 0) return NULL;
    double **A = malloc(m * sizeof(double *));
    if (A == NULL) return NULL;
    A[0] = malloc(m*n*sizeof(double));
    if (A[0] == NULL) {
        free(A);
        return NULL;
 }
    for (int i = 1; i < m; i++)
        A[i] = A[0] + i * n;
    return A;}
// Free it again
void free_2d(double **A) {
 free(A[0]);
 free(A);
};


void main()
{
int i,j,n,m,k;
m = 3;
n = 2;
k = 5;
double **A = malloc_2d(m,k);
for(i = 0; i<m; i++){
    for(j = 0; j<k; j++){
        A[i][j] = 10.0*(i+1)+(j+1);
    }
}

double **B = malloc_2d(k,n);
for(i = 0; i<k; i++){
    for(j = 0; j<n; j++){
        B[i][j] = 20.0*(i+1)+(j+1);
    }
}

double **C = malloc_2d(m,n);

matmult_nat(m,n,k,A,B,C);
// matmult_mkn(m,n,k,m1,m2,m3);
// matmult_nkm(m,n,k,m1,m2,m3);
// matmult_nmk(m,n,k,m1,m2,m3);
// matmult_knm(m,n,k,m1,m2,m3);
// matmult_kmn(m,n,k,m1,m2,m3);
for (int i=0;i<m;i++){
    for (int p=0;p<n;p++){
        printf("%f",C[i][p]);
    }
}


// Deallocating
free_2d(A);
free_2d(B);
free_2d(C);
}