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

    int j = blockIdx.x*blockDim.x+threadIdx.x; //looping through n
    int i = blockIdx.y*blockDim.y+threadIdx.y; //looping through m
    double tmp;
    
    if((j<n)&&(i<m)){
        for(int p=0; p<k; p++){
            //read row of A and col of B 
            //row of A is A[mit*k+kit]
            //col of B is B[kit*n+nit]
            tmp += A[i*k+p] * B[p*n+j];
        }
        //C is C[mit*n+nit]
        C[i*n+j] = tmp;
    }
}

__global__ void
matmult_kernel3(int m, int n, int k, double *A, double *B, double *C){
    //compute C(i,j) and C(i,j+1)
    int j = 2*(blockIdx.x*blockDim.x+threadIdx.x); //looping through n (only half as many threads/blocks)
    int i = blockIdx.y*blockDim.y+threadIdx.y;     //looping through m
    double tmp[2];

    if((j<n)&&(i<m)){

        //additional j to compute (here, either 1 or 0)
        int j_add = MIN(1,n-1-j);

        #pragma unroll
        for(int u=0; u<=j_add; u++){
            tmp[u] = 0.0;
        }
        for(int p=0; p<k; p++){
            //row of A and col of B
            #pragma unroll
            for(int u=0; u<=j_add; u++){
                tmp[u] += A[i*k+p] * B[p*n+j+u];
            }
        }
        //#pragma unroll
        for(int u=0;  u<=j_add; u++){
            C[i*n+j+u] = tmp[u];
        }
        
    }
    
}

#define MAX_ELEM_PER_THREAD 32
__global__ void
matmult_kernel4(int m, int n, int k, double *A, double *B, double *C, int s){
    //compute C(i,j), C(i,j+1), ... C(i,j+s)

    int j = s*(blockIdx.x*blockDim.x+threadIdx.x); //looping through n  (only 1/s as many threads/blocks)
    int i = blockIdx.y*blockDim.y+threadIdx.y; //looping through m
    double tmp[MAX_ELEM_PER_THREAD];

    if((j<n)&&(i<m)){
        //additional j to compute (here, from 0 to s-1)
        int j_add = MAX_ELEM_PER_THREAD-1;//MIN(s-1,n-1-j);
        int real_j_add = MIN(s-1,n-1-j);
        #pragma unroll
        for(int u=0; u<=j_add; u++){
            if(u<real_j_add)
                tmp[u] = 0.0;
        }

        for(int p=0; p<k; p++){
            //row of A and col of B
            #pragma unroll
            for(int u=0; u<=j_add; u++){
                if(u<real_j_add)
                    tmp[u] += A[i*k+p] * B[p*n+j+u];
            }
        }
        #pragma unroll
        for(int u=0;  u<=j_add; u++){
            if(u<real_j_add)
                C[i*n+j+u] = tmp[u];
        }
    }
}




/*
BEGIN GPU 5
##############################################################################
*/

#define BLOCK_SIZE 16

typedef struct {
    int width;
    int height;
    int stride; 
    double* elements;
} Matrix;


// Get a matrix element
__device__ double GetElement(const Matrix A, int row, int col)
{
    return A.elements[row * A.stride + col];
}

// Set a matrix element
__device__ void SetElement(Matrix A, int row, int col,
                           double value)
{
    A.elements[row * A.stride + col] = value;
}

// Get the BLOCK_SIZExBLOCK_SIZE sub-matrix Asub of A that is
// located col sub-matrices to the right and row sub-matrices down
// from the upper-left corner of A
 __device__ Matrix GetSubMatrix(Matrix A, int row, int col) 
{
    Matrix Asub;
    Asub.width    = BLOCK_SIZE;
    Asub.height   = BLOCK_SIZE;
    Asub.stride   = A.stride;
    Asub.elements = &A.elements[A.stride * BLOCK_SIZE * row + BLOCK_SIZE * col];
    return Asub;
}

// Matrix multiplication kernel called by MatMul()
__global__ void gpu5_kernel(const Matrix A, const Matrix B, Matrix C)
{
    // Block row and column
    int blockRow = blockIdx.y;
    int blockCol = blockIdx.x;

    // Each thread block computes one sub-matrix Csub of C
    Matrix Csub = GetSubMatrix(C, blockRow, blockCol);

    // Each thread computes one element of Csub
    // by accumulating results into Cvalue
    double Cvalue = 0;

    // Thread row and column within Csub
    int row = threadIdx.y;
    int col = threadIdx.x;

    // Loop over all the sub-matrices of A and B that are
    // required to compute Csub
    // Multiply each pair of sub-matrices together
    // and accumulate the results
    for (int m = 0; m < (A.width / BLOCK_SIZE); ++m) {

        // Get sub-matrix Asub of A
        Matrix Asub = GetSubMatrix(A, blockRow, m);

        // Get sub-matrix Bsub of B
        Matrix Bsub = GetSubMatrix(B, m, blockCol);

        // Shared memory used to store Asub and Bsub respectively
        __shared__ double As[BLOCK_SIZE][BLOCK_SIZE];
        __shared__ double Bs[BLOCK_SIZE][BLOCK_SIZE];

        // Load Asub and Bsub from device memory to shared memory
        // Each thread loads one element of each sub-matrix
        As[row][col] = GetElement(Asub, row, col);
        Bs[row][col] = GetElement(Bsub, row, col);

        // Synchronize to make sure the sub-matrices are loaded
        // before starting the computation
        __syncthreads();
        // Multiply Asub and Bsub together
        for (int e = 0; e < BLOCK_SIZE; ++e)
            Cvalue += As[row][e] * Bs[e][col];

        // Synchronize to make sure that the preceding
        // computation is done before loading two new
        // sub-matrices of A and B in the next iteration
        __syncthreads();
    }

    // Write Csub to device memory
    // Each thread writes one element
    SetElement(Csub, row, col, Cvalue);
}