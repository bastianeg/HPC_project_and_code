/* jacobi.c - Poisson problem in 3d
 *
 */
#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>
#include <stdio.h>

__global__ void 
initmat(int jmp, double* U, double* Uold, double* F){
    double deltasq = 4.0/((double) (jmp-2) * (double) (jmp-2));
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    if(i<jmp*jmp*jmp){
        Uold[i] = U[i];
        F[i] *= deltasq;
    }
}

__global__ void 
updmat(int jmp, double* U, double* Uold){

    int i = blockIdx.x*blockDim.x+threadIdx.x;
    if(i<jmp*jmp*jmp){
        Uold[i] = U[i];
    }
}

__global__ void 
jacgpu(int jmp, double* U, double* Uold,double* F){
    double onesixth = 1.0/6.0;
    int i = blockIdx.x*blockDim.x+threadIdx.x+1;
    int j = blockIdx.y*blockDim.y+threadIdx.y+1;
    int k = blockIdx.z*blockDim.z+threadIdx.z+1;
    int idx=i+j*jmp+k*jmp*jmp;
    U[idx] = onesixth*(Uold[idx-1]+Uold[idx+1]+Uold[idx-jmp]+\
    Uold[idx+jmp]+Uold[idx+jmp*jmp]+Uold[idx-jmp*jmp]+F[idx]);

}

void
jacobinaive(double **U, double *F, double **Uold, int N, int iter_max,double **tmp) {
    int B=10; // Block size
    double ts, te; // for timing
    double *tmp;
    //define norm and max_iter and Uold and iter and threshold
    int iter = 0;
    int jmp = N+2;

    // update Uold = U
    initmat<<<jmp*jmp*jmp/B,B>>>(jmp, *U, *Uold,F);
    cudaDeviceSynchronize();

    ts = omp_get_wtime();
    //while condition is not satisfied
    while(iter < iter_max)
    {
        jacgpu<<<dim3(N/B,N/B,N/B),dim3(B,B,B)>>>(jmp, *U, *Uold,F);
        cudaDeviceSynchronize();

        // update iteration and Uold
        iter ++;
        *tmp = *U;
        *U = *Uold;
        *Uols = *tmp;

        //updmat<<<jmp*jmp*jmp/B,B>>>(jmp, U,Uold);
        //cudaDeviceSynchronize();
    }
    te = omp_get_wtime() - ts;
    
    printf("Number of iterations: %d\n", iter);
    printf("Elapsed time: %lf\n", te);
    printf("Iterations per second: %lf\n", iter/te);
    printf("Gflop/s: %.3f\n",1e-9*7*N*N*N*iter/te);
    //printf("%.5lf, %.5lf\n", te, 1e-6*11*N*N*N*iter/te);
}
