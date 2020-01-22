#include <stdio.h>
__global__ void
matmult_kernel1(int m, int n, int k, double *A, double *B, double *C){
    
    
    for(int i = 0; i<k; i++){
            for(int j = 0; j<n; j++){
                printf("%.2lf ",B[i*n+j]);
            }
            printf("\n");
        }
        printf("is\n");
    
    // set C to zeros
    for (int i=0;i<m;i++){
        for (int p=0;p<n;p++){
            C[i*n+p]=0; //C[i][p]
        }
    }
    // do matmult with mkn loop order
    for (int i=0;i<m;i++) {
        for (int j=0;j<k;j++){
            for (int p=0;p<n;p++) { 
                C[i*n+p] += A[i*k+j] * B[j*n+p];
            }
        }
    }

    for(int i = 0; i<m; i++){
            for(int j = 0; j<n; j++){
                printf("%.2lf ",C[i*n+j]);
            }
            printf("\n");
        }

}