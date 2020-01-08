#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include "matmult_lib.h"


double** malloc_2d(int m, int n){
    int i;
    double **A = malloc(m * sizeof(double *));
    A[0] = malloc(m*n*sizeof(double));
    for(i = 1; i<m; i++){
        A[i] = A[0] + i*n;
    }

    return A;
}

void free_2d(double **A){
    free(A[0]);
    free(A);
}

int main(void){
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
  
    // Instead of calling cblas_dgemm we can now use this wrapper function:
    matmult_lib(m,n,k,A,B,C);
    
    for(i = 0; i<m; i++){
        for(j = 0; j<k; j++){
            printf("%f ",A[i][j]);
        }
        printf("\n");
    }
    printf("times\n");
    for(i = 0; i<k; i++){
        for(j = 0; j<n; j++){
            printf("%f ",B[i][j]);
        }
        printf("\n");
    }
    printf("equals\n");
    for(i = 0; i<m; i++){
        for(j = 0; j<n; j++){
            printf("%f ",C[i][j]);
        }
        printf("\n");
    }
    free_2d(A);
    free_2d(B);
    free_2d(C);
    return 0;
}