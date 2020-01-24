
__global__ void matmult_kernel1(int m, int n, int k, double *A, double *B, double *C);
__global__ void matmult_kernel2(int m, int n, int k, double *A, double *B, double *C);
__global__ void matmult_kernel3(int m, int n, int k, double *A, double *B, double *C);
__global__ void matmult_kernel4(int m, int n, int k, double *A, double *B, double *C, int s);



// Matrices are stored in row-major order:
// M(row, col) = *(M.elements + row * M.stride + col)
typedef struct {
    int width;
    int height;
    int stride; 
    double* elements;
} Matrix;


#define BLOCK_SIZE 16
__device__ double GetElement(const Matrix A, int row, int col);
__device__ void SetElement(Matrix A, int row, int col, double value);
 __device__ Matrix GetSubMatrix(Matrix A, int row, int col);
__global__ void gpu5_kernel(Matrix A, Matrix B, Matrix C);