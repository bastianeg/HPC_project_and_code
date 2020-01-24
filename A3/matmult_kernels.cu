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
    double tmp1 = 0.0;
    double tmp2 = 0.0;

    if((j<n)&&(i<m)){

        //additional j to compute (here, either 1 or 0)
        int j_add = MIN(1,n-1-j);

        for(int p=0; p<k; p++){
            //row of A and col of B
            tmp1 += A[i*k+p] * B[p*n+j];
            if(j_add == 1)
                tmp2 += A[i*k+p] * B[p*n+j+1];
        }

        C[i*n+j] = tmp1;
        if(j_add == 1)
            C[i*n+j+1] = tmp2;

    }
    
}

__global__ void
matmult_kernel4(int m, int n, int k, double *A, double *B, double *C, int s){
    //compute C(i,j), C(i,j+1), ... C(i,j+s)

    int j = s*(blockIdx.x*blockDim.x+threadIdx.x); //looping through n  (only 1/s as many threads/blocks)
    int i = blockIdx.y*blockDim.y+threadIdx.y; //looping through m
    double tmp1 = 0.0;
    double tmp2 = 0.0;
    double tmp3 = 0.0;
    double tmp4 = 0.0;
    double tmp5 = 0.0;
    double tmp6 = 0.0;
    double tmp7 = 0.0;
    double tmp8 = 0.0;
    double tmp9 = 0.0;
    double tmp10 = 0.0;
    double tmp11 = 0.0;
    double tmp12 = 0.0;
    double tmp13 = 0.0;
    double tmp14 = 0.0;
    double tmp15 = 0.0;
    double tmp16 = 0.0;

    if((j<n)&&(i<m)){
        //additional j to compute (here, from 0 to s-1)
        int j_add = MIN(s-1,n-1-j);

        for(int p=0; p<k; p++){
            //row of A and col of B
            tmp1 += A[i*k+p] * B[p*n+j+0];
            if(j_add > 0)
                tmp2 += A[i*k+p] * B[p*n+j+1];
                if(j_add > 1)
                    tmp3 += A[i*k+p] * B[p*n+j+2];
                    if(j_add > 2)
                        tmp4 += A[i*k+p] * B[p*n+j+3];
                        if(j_add > 3)
                            tmp5 += A[i*k+p] * B[p*n+j+4];
                            if(j_add > 4)
                                tmp6 += A[i*k+p] * B[p*n+j+5];
                                if(j_add > 5)
                                    tmp7 += A[i*k+p] * B[p*n+j+6];
                                    if(j_add > 6)
                                        tmp8 += A[i*k+p] * B[p*n+j+7];
                                        if(j_add > 7)
                                            tmp9 += A[i*k+p] * B[p*n+j+8];
                                            if(j_add > 8)
                                                tmp10 += A[i*k+p] * B[p*n+j+9];
                                                if(j_add > 9)
                                                    tmp11 += A[i*k+p] * B[p*n+j+10];
                                                    if(j_add > 10)
                                                        tmp12 += A[i*k+p] * B[p*n+j+11];
                                                        if(j_add > 11)
                                                            tmp13 += A[i*k+p] * B[p*n+j+12];
                                                            if(j_add > 12)
                                                                tmp14 += A[i*k+p] * B[p*n+j+13];
                                                                if(j_add > 13)
                                                                    tmp15 += A[i*k+p] * B[p*n+j+14];
                                                                    if(j_add > 14)
                                                                        tmp16 += A[i*k+p] * B[p*n+j+15];
        }
        C[i*n+j] = tmp1;
            if(j_add > 0)
                C[i*n+j+1] = tmp2;
                if(j_add > 1)
                    C[i*n+j+2] = tmp3;
                    if(j_add > 2)
                        C[i*n+j+3] = tmp4;
                        if(j_add > 3)
                            C[i*n+j+4] = tmp5;
                            if(j_add > 4)
                                C[i*n+j+5] = tmp6;
                                if(j_add > 5)
                                    C[i*n+j+6] = tmp7;
                                    if(j_add > 6)
                                        C[i*n+j+7] = tmp8;
                                        if(j_add > 7)
                                            C[i*n+j+8] = tmp9;
                                            if(j_add > 8)
                                                C[i*n+j+9] = tmp10;
                                                if(j_add > 9)
                                                    C[i*n+j+10] = tmp11;
                                                    if(j_add > 10)
                                                        C[i*n+j+11] = tmp12;
                                                        if(j_add > 11)
                                                            C[i*n+j+12] = tmp13;
                                                            if(j_add > 12)
                                                                C[i*n+j+13] = tmp14;
                                                                if(j_add > 13)
                                                                    C[i*n+j+14] = tmp15;
                                                                    if(j_add > 14)
                                                                        C[i*n+j+15] = tmp16;
        
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