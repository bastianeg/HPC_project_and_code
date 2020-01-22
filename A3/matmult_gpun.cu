#include "matmult_kernels.h"
#include <stdio.h>
extern "C"{

    void matmult_gpu1(int m, int n, int k, double *A, double *B, double *C){
        
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
        matmult_kernel1<<<1,1>>>(m, n, k, d_A, d_B, d_C);
        cudaDeviceSynchronize();

        //move C back to host
        cudaMemcpy(C, d_C, m*n*sizeof(double), cudaMemcpyDeviceToHost);

        cudaFree(d_A);
        cudaFree(d_B);
        cudaFree(d_C);

    }

    void matmult_gpu2(int m, int n, int k, double *A, double *B, double *C){
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

        //number of blocks is ceil of N/bs 
        int bs = 32;
        int mblocks = m/bs + (int) (m%bs!=0);
        int nblocks = n/bs + (int) (n%bs!=0);

        //call kernel
        matmult_kernel2<<<dim3 (mblocks,nblocks),dim3 (bs,bs)>>>(m, n, k, d_A, d_B, d_C);
        cudaDeviceSynchronize();

        //move C back to host
        cudaMemcpy(C, d_C, m*n*sizeof(double), cudaMemcpyDeviceToHost);

        cudaFree(d_A);
        cudaFree(d_B);
        cudaFree(d_C);
    }

    void matmult_gpu3(int m, int n, int k, double *A, double *B, double *C){
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

        //number of blocks is ceil of N/bs 
        int bs = 32;
        int mblocks = m/bs + (int) (m%bs!=0);
        int nblocks = n/bs/2 + (int) (n%(bs*2)!=0);

        //call kernel
        matmult_kernel3<<<dim3 (mblocks,nblocks),dim3 (bs,bs)>>>(m, n, k, d_A, d_B, d_C);
        cudaDeviceSynchronize();

        //move C back to host
        cudaMemcpy(C, d_C, m*n*sizeof(double), cudaMemcpyDeviceToHost);

        cudaFree(d_A);
        cudaFree(d_B);
        cudaFree(d_C);
    }

    void matmult_gpu4(int m, int n, int k, double *A, double *B, double *C){
        //number of elements to compute in each thread
        int s = (int) getenv("NUM_ELEM_PER_THREAD");
        printf("%d\n",s);
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

        //number of blocks is ceil of N/bs 
        int bs = 32;
        int mblocks = m/bs + (int) (m%bs!=0);
        int nblocks = n/bs/s + (int) (n%(bs*s)!=0);

        //call kernel
        matmult_kernel4<<<dim3 (mblocks,nblocks),dim3 (bs,bs)>>>(m, n, k, d_A, d_B, d_C, s);
        cudaDeviceSynchronize();

        //move C back to host
        cudaMemcpy(C, d_C, m*n*sizeof(double), cudaMemcpyDeviceToHost);

        cudaFree(d_A);
        cudaFree(d_B);
        cudaFree(d_C);
    }

}