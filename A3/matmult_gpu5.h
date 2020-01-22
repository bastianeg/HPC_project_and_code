
__device__ float GetElement(const Matrix A, int row, int col);
__device__ void SetElement(Matrix A, int row, int col, float value);
 __device__ Matrix GetSubMatrix(Matrix A, int row, int col);
__global__ void MatMulKernel(double* A, double* B, double* C);
void MatMul(int m, int n, int k, const Matrix A, const Matrix B, Matrix C);


