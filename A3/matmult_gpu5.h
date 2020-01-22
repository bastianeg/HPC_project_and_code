typedef struct {
    int width;
    int height;
    int stride; 
    float* elements;
} Matrix;
__device__ float GetElement(const Matrix A, int row, int col);
__device__ void SetElement(Matrix A, int row, int col, float value);
 __device__ Matrix GetSubMatrix(Matrix A, int row, int col);
__global__ void MatMulKernel5(double* A, double* B, double* C);


