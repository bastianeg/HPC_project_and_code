#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

void matmult_lib(int m, int n, int k, double **A, double **B, double **C){

    cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans,m,n,k,1.0,*A,k,*B,n,0.0,*C,n);
}