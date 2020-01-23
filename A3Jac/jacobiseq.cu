/* jacobi.c - Poisson problem in 3d
 *
 */
#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>
#include <stdio.h>

__global__ void 
initmat(int N, double* U, double* Uold, double* F,double deltasq){
    int jmp=N+2;
    int idx;
    for(int i = 1; i<(N+1); i++){
        for(int j = 1; j<(N+1); j++){
            for(int k = 1; k<(N+1); k++){
                idx = jmp*jmp*i+jmp*j+k;
                Uold[idx] = U[idx];
                F[idx] *= deltasq;
            }
        }
    }
    
}

__global__ void 
updmat(int N, double* U, double* Uold){
    int jmp=N+2;
    int idx;
    for(int i = 1; i<(N+1); i++){
        for(int j = 1; j<(N+1); j++){
            for(int k = 1; k<(N+1); k++){
                idx = jmp*jmp*i+jmp*j+k;
                Uold[idx] = U[idx];
            }
        }
    }
}

__global__ void 
jacgpu(int N, double* U, double* Uold, double* F){
    int jmp=N+2;
    int idx;
    double onesixth = 1.0/6.0;
    for(int i = 1; i<(N+1); i++){
        for(int j = 1; j<(N+1); j++){
            for(int k = 1; k<(N+1); k++){
                idx = jmp*jmp*i+jmp*j+k;
                U[idx] = onesixth*(Uold[idx-1]+Uold[idx+1]+Uold[idx-jmp]\
                +Uold[idx+jmp]+Uold[idx-jmp*jmp]+Uold[idx+jmp*jmp]); //+F[idx]
            }
        }
    }
}

void
jacobiseq(double *U, double *F, double *Uold, int N, int iter_max) {
    double ts, te; // for timing
    double deltasq = 4.0/((double) N * (double) N);
    //define norm and max_iter and Uold and iter and threshold
    int iter = 0;

    // update Uold = U
    initmat<<<1,1>>>(N, U,Uold,F,deltasq);
    cudaDeviceSynchronize();

    ts = omp_get_wtime();
    //while condition is not satisfied
    while((iter < iter_max))
    {
        jacgpu<<<1,1>>>(N, U, Uold, F);
        cudaDeviceSynchronize();

        // update iteration and Uold
        iter ++;

        //pointer swap instead
        updmat<<<1,1>>>(N, U, Uold);
        cudaDeviceSynchronize();

    }
    te = omp_get_wtime() - ts;
    
    printf("Number of iterations: %d\n", iter);
    printf("Elapsed time: %lf\n", te);
    printf("Iterations per second: %lf\n", iter/te);
    //printf("%.5lf, %.5lf\n", te, 1e-6*11*N*N*N*iter/te);
}
