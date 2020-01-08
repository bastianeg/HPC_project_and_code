#include <Accelerate/Accelerate.h>

enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
enum CBLAS_UPLO {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE {CblasLeft=141, CblasRight=142};

enum CBLAS_ORDER order = CblasColMajor;
enum CBLAS_TRANSPOSE transA = CblasNoTrans;
enum CBLAS_TRANSPOSE transB = CblasNoTrans;


void cblas_dgemm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TA,
                 const enum CBLAS_TRANSPOSE TB,
                 const int M, const int N, const int K,
                 const double  alpha, const double *A, const int lda,
                 const double *B, const int ldb, const double  beta,
                 double *C, const int ldc);

void matmult_lib(int m, int n, int k, double **A, double **B, double **C){
       
    // cblas_dgemm(order,transA,transB, rowsA, colsB, common ,1.0, A, rowsA ,B, common ,0.0,C, rowsA);
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans,m,n,k,1.0,*A,k,*B,n,0.0,*C,n);
}