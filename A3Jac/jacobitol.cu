/* jacobi.c - Poisson problem in 3d
 *
 */
 #ifdef _OPENMP
 #include <omp.h>
 #endif
 
 #include <math.h>
 #include <stdio.h>
 

 __inline__ __device__
 double warpReduceSum(double value) {
    //printf("%f  ",value);
    for (int i = 16; i > 0; i /= 2){
        //printf("%i  ",i);
        //printf("%f  ",value);
        value += __shfl_down_sync(-1, value, i);
    }
    return value;
 }

 __inline__ __device__
double blockReduceSum(double value) {
    __shared__ double smem[32]; // Max 32 warp sums
    if (threadIdx.x < warpSize)
        smem[threadIdx.x] = 0;

    __syncthreads();
    value = warpReduceSum(value);
    
    if (threadIdx.x % warpSize == 0)
        smem[threadIdx.x / warpSize] = value;

    __syncthreads();
    
    if (threadIdx.x < warpSize)
        value = smem[threadIdx.x];
    
    //printf("Reenter %d",threadIdx.x);
    return warpReduceSum(value);
}

 __global__ void 
 reduction_presum (double *dpart, int n, double* res)
 {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    double value = 0.0;
    for (int i = idx; i < n; i += blockDim.x * gridDim.x){
        value += dpart[i];
    }
    value = idx < n ? value : 0;
    value = blockReduceSum(value);
    if (threadIdx.x == 0){
         atomicAdd(res, value);
    }
 }
 /*
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
*/
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
 diff(int jmp,double* dpart){ // double* U, double* Uold
 
     int i = blockIdx.x*blockDim.x+threadIdx.x;
     if(i<jmp*jmp*jmp){
         dpart[i] =1; // U[i]-Uold[i];
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
 jacobitol(double *U, double *F, double *Uold, int N, int iter_max, double tol,double* dpart) { 
     int B=10; // Block size
     double ts, te; // for timing
     //define norm and max_iter and Uold and iter and threshold
     int iter = 0;
     int jmp = N+2;
     double d=tol+8.0;
     double res=0;
     // update Uold = U
     initmat<<<jmp*jmp*jmp/B,B>>>(jmp, U,Uold,F);
     cudaDeviceSynchronize();

     ts = omp_get_wtime();
     //while condition is not satisfied
     while(iter<iter_max) //(d>tol) && (iter < iter_max))
     {
         res = 0.0;
         jacgpu<<<dim3(N/B,N/B,N/B),dim3(B,B,B)>>>(jmp, U, Uold,F);
         cudaDeviceSynchronize();
         
         //Calculate d
         diff<<<jmp*jmp*jmp/B,B>>>(jmp,dpart);
         cudaDeviceSynchronize();

         reduction_presum<<<jmp*jmp*jmp/B,B>>>(dpart, jmp*jmp*jmp, &res);
         cudaDeviceSynchronize();
         printf("d: %f\n",res);
         //printf("%f",res);
         //update iteration and Uold
         iter ++;

         updmat<<<jmp*jmp*jmp/B,B>>>(jmp, U,Uold);
         cudaDeviceSynchronize();
     }
     te = omp_get_wtime() - ts;
     
     printf("Number of iterations: %d\n", iter);
     printf("Elapsed time: %lf\n", te);
     printf("Iterations per second: %lf\n", iter/te);
     printf("Gflop/s: %.3f\n",1e-9*7*N*N*N*iter/te);
     //printf("%.5lf, %.5lf\n", te, 1e-6*11*N*N*N*iter/te);
 }
 