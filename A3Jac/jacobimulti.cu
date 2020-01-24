/* jacobi.c - Poisson problem in 3d
 *
 */
#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>
#include <stdio.h>

__global__ void 
initmat(int jmp, double* U, double* Uold, double* F, double deltasq){

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
jaclower(int jmp, double* U, double* Uold, double* upper_Uold, double* F, double onesixth){

    int i = blockIdx.x*blockDim.x+threadIdx.x+1; // goes from 1 to N
    int j = blockIdx.y*blockDim.y+threadIdx.y+1; // goes from 1 to N
    int k = blockIdx.z*blockDim.z+threadIdx.z+1; // goes from 1 to N/2  0=1 

    int idx=i + j*jmp + k*jmp*jmp;
    if((i<(jmp-1)) && (j<(jmp-1)) && (k<((jmp-2)/2)+1)){
    if( k == ((jmp-2)/2) ){
        U[idx] = onesixth*(Uold[idx-1]+Uold[idx+1]+Uold[idx-jmp]+\
        Uold[idx+jmp]+Uold[idx-jmp*jmp]+upper_Uold[i+j*jmp]+F[idx]);
    }else{
        U[idx] = onesixth*(Uold[idx-1]+Uold[idx+1]+Uold[idx-jmp]+\
        Uold[idx+jmp]+Uold[idx+jmp*jmp]+Uold[idx-jmp*jmp]+F[idx]);
    }
    }
}

__global__ void 
jacupper(int jmp, double* U, double* Uold,double* lower_Uold, double* F, double onesixth){

    int i = blockIdx.x*blockDim.x+threadIdx.x+1; // goes from 1 to N
    int j = blockIdx.y*blockDim.y+threadIdx.y+1; // goes from 1 to N
    int k = blockIdx.z*blockDim.z+threadIdx.z;   // goes from 0 to N/2-1   9=4
    if((i<(jmp-1)) && (j<(jmp-1)) && (k<(((jmp-2)/2)))){
    int idx=i + j*jmp + k*jmp*jmp;

    if(k==0){
        U[idx] = onesixth*(Uold[idx-1]+Uold[idx+1]+Uold[idx-jmp]+\
        Uold[idx+jmp]+Uold[idx+jmp*jmp]+lower_Uold[i+j*jmp+((jmp-2)/2)*jmp*jmp]+F[idx]);

    }else{
        U[idx] = onesixth*(Uold[idx-1]+Uold[idx+1]+Uold[idx-jmp]+\
        Uold[idx+jmp]+Uold[idx+jmp*jmp]+Uold[idx-jmp*jmp]+F[idx]);
    }

}
}

void
jacobimulti(double* D0U,double* D1U, double* D0F, double* D1F, double* D0Uold, double* D1Uold, int N, int iter_max) {
    int B=10; // Block size 2B or not 2B?

    double ts, te; // for timing
    double deltasq = 4.0/((double) N * (double) N);
    //define norm and max_iter and Uold and iter and threshold
    int iter = 0;
    double onesixth = 1.0/6.0;
    int jmp = N + 2;
<<<<<<< HEAD
    int halfjmp=jmp/2 + (int) (N%B!=0);
    int nhalf_blocks=N/(2*B) + (int) (N%(2*B)!=0);
=======
    int halfjmp=jmp/(B*2) + (int) (jmp%(2*B)!=0);
    int halfN=N/(2*B) + (int) (N%(2*B)!=0);
>>>>>>> a891502022529a65948c3b5281daca06f114281e
    int n_blocks = N/B + (int) (N%B!=0);
    int jmp_blocks = jmp*jmp*halfjmp/(B*B*B) + (int) (jmp%B!=0);
    // update Uold = U
    cudaSetDevice(0);
    initmat<<<jmp_blocks,(B*B*B)>>>(jmp, D0U, D0Uold, D0F, deltasq);

    cudaSetDevice(1);
    initmat<<<jmp_blocks,B>>>(jmp, D1U, D1Uold, D1F, deltasq);
    cudaDeviceSynchronize();

    ts = omp_get_wtime();
    //while condition is not satisfied
    while((iter < iter_max))
    {

        cudaSetDevice(0);
        jaclower<<<dim3(n_blocks,n_blocks,nhalf_blocks),dim3(B,B,B)>>>(jmp, D0U, D0Uold, D1Uold ,D0F, onesixth);
        cudaSetDevice(1);
        jacupper<<<dim3(n_blocks,n_blocks,nhalf_blocks),dim3(B,B,B)>>>(jmp, D1U, D1Uold, D0Uold, D1F, onesixth);
        cudaDeviceSynchronize();

        // update iteration and Uold
        iter ++;

        cudaSetDevice(0);
        updmat<<<jmp_blocks,(B*B*B)>>>(jmp, D0U, D0Uold);

        cudaSetDevice(1);
        updmat<<<jmp_blocks,(B*B*B)>>>(jmp, D1U, D1Uold);
        cudaDeviceSynchronize();
}
    te = omp_get_wtime() - ts;
    
    printf("Number of iterations: %d\n", iter);
    printf("Elapsed time: %lf\n", te);
    printf("Iterations per second: %lf\n", iter/te);
    printf("Gflop/s: %.3f\n",1e-9*7*N*N*N*iter/te);
    //printf("%.5lf, %.5lf\n", te, 1e-6*11*N*N*N*iter/te);
}
