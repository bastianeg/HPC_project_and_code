/* jacobi.c - Poisson problem in 3d
 *
 */
#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>
#include <stdio.h>

__global__ void 
initmat(int N, double* U, double* Uold, double* F, double deltasq){

    int i = blockIdx.x*blockDim.x+threadIdx.x;
    Uold[i] = U[i];
    F[i] *= deltasq;

}

__global__ void 
updmat(int N, double* U, double* Uold){

    int i = blockIdx.x*blockDim.x+threadIdx.x;
    Uold[i] = U[i];

}

__global__ void 
jacgpu(int N, double* U, double* Uold,double* F, double onesixth){

    int i = blockIdx.x*blockDim.x+threadIdx.x;
    int j = blockIdx.y*blockDim.y+threadIdx.y;
    int k = blockIdx.z*blockDim.z+threadIdx.z;
    int jmp=N+2;
    int idx=i+j*jmp+k*jmp*jmp;
    U[idx] = onesixth*(Uold[idx-1]+Uold[idx+1]+Uold[idx-jmp]+\
    Uold[idx+jmp]+Uold[idx+jmp*jmp]+Uold[idx-jmp*jmp]+F[idx]);

}

void
jacobinaive(double *U, double *F, double *Uold, int N, int iter_max) {
    int B=1; // Block size

    double ts, te; // for timing
    double deltasq = 4.0/((double) N * (double) N);
    //define norm and max_iter and Uold and iter and threshold
    int iter = 0;
    double onesixth = 1.0/6.0;

    // update Uold = U
    initmat<<<N*N*N/B,B>>>(N, U,Uold,F,deltasq);
    cudaDeviceSynchronize();

    ts = omp_get_wtime();
    //while condition is not satisfied
    while(iter < iter_max)
    {
        jacgpu<<<dim3(N/B,N/B,N/B),dim3(B,B,B)>>>(N, U, Uold,F, onesixth);
        cudaDeviceSynchronize();

        // update iteration and Uold
        iter ++;

        updmat<<<N*N*N/B,B>>>(N, U,Uold);
        cudaDeviceSynchronize();
    }
    te = omp_get_wtime() - ts;
    
    printf("Number of iterations: %d\n", iter);
    printf("Elapsed time: %lf\n", te);
    printf("Iterations per second: %lf\n", iter/te);
    //printf("%.5lf, %.5lf\n", te, 1e-6*11*N*N*N*iter/te);
}
