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

    int i = threadIdx.x;
    Uold[i] = U[i];
    F[i] *= deltasq;

}

__global__ void 
updmat(int N, double* U, double* Uold){

    int i = threadIdx.x;
    Uold[i] = U[i];

}

__global__ void 
jacgpupper(int N, double* U, double* Uold,double* F, double onesixth){

    int i = threadIdx.x;
    int j = threadIdx.y;
    int k = threadIdx.z;
    int jmp=N+2;
    int Nj=jmp*j;
    int N2k=jmp*jmp*k;
    int NoWall=jmp*jmp+1;
    U[i+Nj+N2k+NoWall] = onesixth*(Uold[i+Nj+N2k-1+NoWall]+Uold[i+Nj+N2k+1+NoWall]+Uold[i+Nj+N2k-N+NoWall]+\
    Uold[i+Nj+N2k+N+NoWall]+Uold[i+Nj+N2k-N*N+NoWall]+Uold[i+Nj+N2k+N*N+NoWall]+F[i+Nj+N2k+NoWall]);

}

__global__ void 
jacglower(int N, double* U, double* Uold,double* F, double onesixth){

    int i = threadIdx.x;
    int j = threadIdx.y;
    int k = threadIdx.z;
    int jmp=N+2;
    int Nj=jmp*j;
    int N2k=jmp*jmp*k;
    int NoWall=1;
    U[i+Nj+N2k+NoWall] = onesixth*(Uold[i+Nj+N2k-1+NoWall]+Uold[i+Nj+N2k+1+NoWall]+Uold[i+Nj+N2k-N+NoWall]+\
    Uold[i+Nj+N2k+N+NoWall]+Uold[i+Nj+N2k-N*N+NoWall]+Uold[i+Nj+N2k+N*N+NoWall]+F[i+Nj+N2k+NoWall]);

}

void
jacobimulti(double* D0U,double* D1U, double* D0F, double* D1F, double* D0Uold, double* D1old, int N, int iter_max) {
    int B=1; // Block size

    double ts, te; // for timing
    double deltasq = 4.0/((double) N * (double) N);
    //define norm and max_iter and Uold and iter and threshold
    int iter = 0;
    double onesixth = 1.0/6.0;
    int halfN=N/(2*B);

    // update Uold = U
    cudaSetDevice(0);
    initmat<<<N*N*halfN,B>>>(N, D0U, D0Uold, D0F, deltasq);

    cudaSetDevice(1);
    initmat<<<N*N*halfN,B>>>(N, D1U, D1old, D1F, deltasq);
    cudaDeviceSynchronize();

    ts = omp_get_wtime();
    //while condition is not satisfied
    while((iter < iter_max))
    {
        // start wallclock timer

        // from  i to j to k
        // for i
        cudaSetDevice(0);
        jacgpupper<<<dim3(halfN,halfN,halfN),dim3(B,B,B)>>>(N, D0U, D0Uold, D0F, onesixth);
        cudaSetDevice(1);
        jaclower<<<dim3(halfN,halfN,halfN),dim3(B,B,B)>>>(N, D1U, D1Uold, D1F, onesixth);
        cudaDeviceSynchronize();

        // update iteration and Uold
        iter ++;

        cudaSetDevice(0);
        updmat<<<N*N*halfN,B>>>(N, D0U, D0Uold);

        cudaSetDevice(1);
        updmat<<<N*N*N/(2*B),B>>>(N, D1U, D1old);
        cudaDeviceSynchronize();
}
    te = omp_get_wtime() - ts;
    
    printf("Number of iterations: %d\n", iter);
    printf("Elapsed time: %lf\n", te);
    printf("Iterations per second: %lf\n", iter/te);
    //printf("%.5lf, %.5lf\n", te, 1e-6*11*N*N*N*iter/te);
}
