#include "matmult_kernels.h"
#include <stdio.h>
extern "C"{

    void matmult_gpu1(int m, int n, int k, double *A, double *B, double *C){
        for(int i = 0; i<m; i++){
            for(int j = 0; j<k; j++){
                printf("%.2lf ",A[i*n+j]);
            }
            printf("\n");
        }
        printf("times\n");
        for(int i = 0; i<k; i++){
            for(int j = 0; j<n; j++){
                printf("%.2lf ",B[i*n+j]);
            }
            printf("\n");
        }
        printf("is\n");
        
        //allocate memory on GPU
        double* d_A;
        double* d_B;
        double* d_C;
        cudaMalloc((void**) &d_A, m*k*sizeof(double));
        cudaMalloc((void**) &d_B, n*k*sizeof(double));
        cudaMalloc((void**) &d_C, m*n*sizeof(double));

        //move A and B to GPU
        cudaMemcpy(d_A, A, m*k*sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_B, B, n*k*sizeof(double), cudaMemcpyHostToDevice);
        
        //call kernel
        matmult_kernel1<<<1,1>>>(m, n, k, A, B, C);
        cudaDeviceSynchronize();

        //move C back to host
        cudaMemcpy(C, d_C, m*n*sizeof(double), cudaMemcpyDeviceToHost);

        cudaFree(d_A);
        cudaFree(d_B);
        cudaFree(d_C);

        for(int i = 0; i<m; i++){
            for(int j = 0; j<n; j++){
                printf("%.2lf ",C[i*n+j]);
            }
            printf("\n");
        }
        printf("on the CPU\n");

    }

    void matmult_gpu2(int m, int n, int k, double *A, double *B, double *C){
        
    }

    void matmult_gpu3(int m, int n, int k, double *A, double *B, double *C){
        
    }

    void matmult_gpu4(int m, int n, int k, double *A, double *B, double *C){
        
    }

}