#include <stdio.h>

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

__global__ void
matmult_kernel1(int m, int n, int k, double *A, double *B, double *C){
    
    // set C to zeros
    for (int i=0;i<m;i++){
        for (int p=0;p<n;p++){
            C[i*n+p]=0; //C[i][p]
        }
    }
    // do matmult with mkn loop order
    for (int i=0;i<m;i++) {
        for (int j=0;j<k;j++){
            for (int p=0;p<n;p++) { 
                C[i*n+p] += A[i*k+j] * B[j*n+p];
            }
        }
    }
}

__global__ void
matmult_kernel2(int m, int n, int k, double *A, double *B, double *C){
    
    int i = blockIdx.x*blockDim.x+threadIdx.x; //looping through m
    int j = blockIdx.y*blockDim.y+threadIdx.y; //looping through n

    if((i<m)&&(j<n)){
        //init C to zero
        C[i*n+j] = 0.0;
        for(int p=0; p<k; p++){
            //read row of A and col of B
            C[i*n+j] += A[i*k+p] * B[p*n+j];
        }
    }
}

__global__ void
matmult_kernel3(int m, int n, int k, double *A, double *B, double *C){
    //compute C(i,j) and C(i,j+1)
    int i = blockIdx.x*blockDim.x+threadIdx.x; //looping through m
    int j = 2*(blockIdx.y*blockDim.y+threadIdx.y); //looping through n (only half as many threads/blocks)
    
    

    if((i<m)&&(j<n)){
        //additional j to compute (here, either 1 or 0)
        int j_add = MIN(1,n-1-j);
        
        //init C to zero
        for(int u=0; u<=j_add; u++){
            C[i*n+j+u] = 0.0;
        }
        C[i*n+j+1] = 0.0;
        for(int p=0; p<k; p++){
            //row of A and col of B
            for(int u=0; u<=j_add; u++){
                C[i*n+j+u] += A[i*k+p] * B[p*n+j+u];
            }
        }
    }
}

__global__ void
matmult_kernel4(int m, int n, int k, double *A, double *B, double *C,int s){
    //compute C(i,j), C(i,j+1), ... C(i,j+s)
    int i = blockIdx.x*blockDim.x+threadIdx.x; //looping through m
    int j = s*(blockIdx.y*blockDim.y+threadIdx.y); //looping through n (only 1/s as many threads/blocks)

    if((i<m)&&(j<n)){
        //additional j to compute (here, from 0 to s-1)
        int j_add = MIN(s-1,n-1-j);
        
        //init C to zero
        for(int u=0; u<=j_add; u++){
            C[i*n+j+u] = 0.0;
        }
        C[i*n+j+1] = 0.0;
        for(int u=0; u<=j_add; u++){
            for(int p=0; p<k; p++){
                //row of A and col of B
                C[i*n+j+u] += A[i*k+p] * B[p*n+j+u];
            }
        }
    }
}