#include "matmult_kernels.h"
#include <stdio.h>
#include "cublas_v2.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>



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
        for(int i=0; i<m; i++){
            for(int j=0; j<k; j++){
                printf("%5.1f ",A[i*k+j]);
            }
            printf("\n");
        }
        printf("times\n");
        for(int i=0; i<k; i++){
            for(int j=0; j<n; j++){
                printf("%5.1f ",B[i*k+j]);
            }
            printf("\n");
        }
        printf("equals\n");
        
        //allocate memory on GPU
        double* d_A;
        double* d_B;
        double* d_C;
        int bs = 16;
        cudaMalloc((void**) &d_A, m*k*sizeof(double));
        cudaMalloc((void**) &d_B, n*k*sizeof(double));
        cudaMalloc((void**) &d_C, m*n*sizeof(double));

        //move A and B to GPU
        cudaMemcpy(d_A, A, m*k*sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_B, B, n*k*sizeof(double), cudaMemcpyHostToDevice);

        //number of blocks is ceil of N/bs
        if(getenv("BLOCK_SIZE")!=NULL){
            bs = atoi(getenv("BLOCK_SIZE"));
        }
        int mblocks = m/bs + (int) (m%bs!=0);
        int nblocks = n/bs + (int) (n%bs!=0);
        
        //call kernel
        matmult_kernel2<<<dim3 (nblocks,mblocks),dim3 (bs,bs)>>>(m, n, k, d_A, d_B, d_C);
        checkCudaErrors(cudaDeviceSynchronize());

        //move C back to host
        cudaMemcpy(C, d_C, m*n*sizeof(double), cudaMemcpyDeviceToHost);

        cudaFree(d_A);
        cudaFree(d_B);
        cudaFree(d_C);

        for(int i=0; i<m; i++){
            for(int j=0; j<n; j++){
                printf("%5.1f ",C[i*k+j]);
            }
            printf("\n");
        }
    }

    void matmult_gpu3(int m, int n, int k, double *A, double *B, double *C){
        //allocate memory on GPU
        double* d_A;
        double* d_B;
        double* d_C;
        int bs = 16;
        cudaMalloc((void**) &d_A, m*k*sizeof(double));
        cudaMalloc((void**) &d_B, n*k*sizeof(double));
        cudaMalloc((void**) &d_C, m*n*sizeof(double));

        //move A and B to GPU
        cudaMemcpy(d_A, A, m*k*sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_B, B, n*k*sizeof(double), cudaMemcpyHostToDevice);

        //number of blocks is ceil of N/bs 
        if(getenv("BLOCK_SIZE")!=NULL){
            int bs = atoi(getenv("BLOCK_SIZE"));
        }
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
        const int s = atoi(getenv("NUM_ELEM_PER_THREAD"));
        //allocate memory on GPU
        double* d_A;
        double* d_B;
        double* d_C;
        int bs = 16;
        cudaMalloc((void**) &d_A, m*k*sizeof(double));
        cudaMalloc((void**) &d_B, n*k*sizeof(double));
        cudaMalloc((void**) &d_C, m*n*sizeof(double));

        //move A and B to GPU
        cudaMemcpy(d_A, A, m*k*sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_B, B, n*k*sizeof(double), cudaMemcpyHostToDevice);

        //number of blocks is ceil of N/bs 
        if(getenv("BLOCK_SIZE")!=NULL){
            bs = atoi(getenv("BLOCK_SIZE"));
        }
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

    void matmult_gpu5(int m, int n, int k, double *A, double *B, double *C){

        // Load A and B to device memory
        Matrix d_A;
        d_A.width = d_A.stride = k; 
        d_A.height = m;
        size_t size = k * m * sizeof(double);
        cudaMalloc((void**) &d_A.elements, size);
        cudaMemcpy(d_A.elements, A, size, cudaMemcpyHostToDevice);

        Matrix d_B;
        d_B.width = d_B.stride = n; 
        d_B.height = k;
        size = n * k * sizeof(double);
        cudaMalloc((void**) &d_B.elements, size);
        cudaMemcpy(d_B.elements, B, size,
        cudaMemcpyHostToDevice);

        // Allocate C in device memory
        Matrix d_C;
        d_C.width = d_C.stride = n; 
        d_C.height = m;
        size = m * n * sizeof(double);
        cudaMalloc((void**) &d_C.elements, size);

        // Invoke kernel
        dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
        dim3 dimGrid(n / dimBlock.x, m / dimBlock.y);
        gpu5_kernel<<<dimGrid, dimBlock>>>(d_A, d_B, d_C);

        // Read C from device memory
        cudaMemcpy(C, d_C.elements, size, cudaMemcpyDeviceToHost);

        // Free device memory
        cudaFree(d_A.elements);
        cudaFree(d_B.elements);
        cudaFree(d_C.elements);
    }


    void matmult_gpulib(int m, int n, int k, double *A, double *B, double *C) {
        
        int lda=m,ldb=k,ldc=m;
        double alf = 1.0;
        double bet = 0.0;
        double *alpha = &alf;
        double *beta = &bet;
        //double *d_alpha;
        //double *d_beta;
        double* d_A;
        double* d_B;
        double* d_C;

        cudaMalloc((void **)&d_A,  m * k * sizeof(double));
        cudaMalloc((void **)&d_B,  k * n * sizeof(double));
        cudaMalloc((void **)&d_C,  n * m * sizeof(double));
        //cudaMalloc((void **)&d_alpha,  sizeof(double));
        //cudaMalloc((void **)&d_beta,  sizeof(double));

        cudaMemcpy(d_A, A, m * k * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_B, B, k * n * sizeof(double), cudaMemcpyHostToDevice);
        //cudaMemcpy(d_alpha, alpha, sizeof(double), cudaMemcpyHostToDevice);
        //cudaMemcpy(d_beta, beta, sizeof(double), cudaMemcpyHostToDevice);

        // Create a handle for CUBLAS
        cublasHandle_t handle;
        cublasCreate(&handle);

        // Do the actual multiplication using the library function
        cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, alpha, d_A, lda, d_B, ldb, beta, d_C, ldc);

        // Destroy the handle
        cublasDestroy(handle);

        // Read C from device memory
        cudaMemcpy(C, d_C, n * m * sizeof(double), cudaMemcpyDeviceToHost);

        // Free device memory
        cudaFree(d_A);
        cudaFree(d_B);
        cudaFree(d_C);
    }

    
}