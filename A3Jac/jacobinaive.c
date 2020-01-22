/* jacobi.c - Poisson problem in 3d
 *
 */
#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>
#include <stdio.h>

__global__ void 
initmat(int n, double* U, double* Uold, double* F,double deltasq){

    int i = threadIdx.x;
    Uold[i] = U[i];
    F[i] *= deltasq;

}

__global__ void 
updmat(int N, double* U, double* Uold){

    int i = threadIdx.x;
    Uold[i] = U[i];
    F[i] *= ;

}

__global__ void 
jacgpu(int N, double* A, double* b, double* result){

    int i = threadIdx.x;
    U[i] = onesixth*(Uold[i-1]+Uold[i+1]+Uold[i+N]+\
    Uold[i-N]+Uold[i+N*N]+Uold[i-N*N]+F[i]);

}

void
jacobi(double *U, double *F, double *Uold, int N, int iter_max, double tol) {
    double ts, te; // for timing
    double deltasq = 4.0/((double) N * (double) N);
    //define norm and max_iter and Uold and iter and threshold
    int iter = 0;
    double onesixth = 1.0/6.0;

    // update Uold = U
    initmat<<<N^2,N>>>(N, U,Uold);
    cudaDeviceSynchronize();

    ts = omp_get_wtime();
    //while condition is not satisfied
    while((iter < iter_max))
    {
        // start wallclock timer

        // from  i to j to k
        // for i
        jacgpu<<<N*N,N>>>(N, U,Uold);
        cudaDeviceSynchronize();

        // update iteration and Uold
        iter ++;

        updmat<<<N*N,N>>>(N, U,Uold);
        cudaDeviceSynchronize();
    }
    te = omp_get_wtime() - ts;
    
    printf("Number of iterations: %d\n", iter);
    printf("Norm: %lf\n", d);
    printf("Elapsed time: %lf\n", te);
    printf("Iterations per second: %lf\n", iter/te);
    //printf("%.5lf, %.5lf\n", te, 1e-6*11*N*N*N*iter/te);
}
