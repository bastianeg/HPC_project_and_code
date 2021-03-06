/* jacobi.c - Poisson problem in 3d
 *
 */
 #ifdef _OPENMP
 #include <omp.h>
 #endif
 
 #include <math.h>
 #include <stdio.h>
 #include <helper_cuda.h>
 

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
 reduction_presum (double *U, double *Uold, int n, double *res)
 {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    double value = 0.0;
    double tmp;
    for (int i = idx; i < n; i += blockDim.x * gridDim.x){
        tmp = U[i]-Uold[i];
        value += tmp*tmp;
    }
    value = idx < n ? value : 0;
    value = blockReduceSum(value);
    if (threadIdx.x == 0){
        atomicAdd(res, value);
    }

 }

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
    if((i<(jmp-1)) && (j<(jmp-1)) && (k<(jmp-1))){
    int idx=i+j*jmp+k*jmp*jmp;
    U[idx] = onesixth*(Uold[idx-1]+Uold[idx+1]+Uold[idx-jmp]+\
    Uold[idx+jmp]+Uold[idx+jmp*jmp]+Uold[idx-jmp*jmp]+F[idx]);
    }
}
 
 void
 jacobitol(double *U, double *F, double *Uold, int N, int iter_max, double tol,double* d_res) { 
     int B=10; // Block size
     double ts, te; // for timing
     //define norm and max_iter and Uold and iter and threshold
     int iter = 0;
     int jmp = N+2;
     double d = N*N*N;
     double res;
     int n_blocks = N/B + (int) (N%B!=0);
     int jmp_blocks = jmp*jmp*jmp/(B*B*B) + (int) (jmp*jmp*jmp%(B*B*B)!=0);

     // update Uold = U
     initmat<<<jmp*jmp*jmp_blocks,B>>>(jmp, U,Uold,F);
     cudaDeviceSynchronize();

     ts = omp_get_wtime();
     //while condition is not satisfied
     while((d/sqrt(N*N*N)>tol) && (iter < iter_max)) //(d>tol) && (iter < iter_max))
     {
         res = 0.0;
         cudaMemcpy(d_res,&res,sizeof(double),cudaMemcpyHostToDevice);
         jacgpu<<<dim3(n_blocks,n_blocks,n_blocks),dim3(B,B,B)>>>(jmp, U, Uold,F);
         cudaDeviceSynchronize();

         reduction_presum<<<jmp_blocks,(B*B*B)>>>(U,Uold, jmp*jmp*jmp, d_res);
         checkCudaErrors(cudaDeviceSynchronize());
         cudaMemcpy(&res,d_res,sizeof(double),cudaMemcpyDeviceToHost);

         d = sqrt(res);
         //update iteration and Uold
         iter ++;

         updmat<<<jmp_blocks,(B*B*B)>>>(jmp, U,Uold);
         cudaDeviceSynchronize();
     }
     te = omp_get_wtime() - ts;
     
     printf("Number of iterations: %d\n", iter);
     printf("Elapsed time: %lf\n", te);
     printf("Iterations per second: %lf\n", iter/te);
     printf("Gflop/s: %.3f\n",1e-9*7*N*N*N*iter/te);
     //printf("%.5lf, %.5lf\n", te, 1e-6*11*N*N*N*iter/te);
 }
 