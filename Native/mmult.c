#include <stdlib.h>
#include <stdio.h>

// allocate a double-prec m x n matrix
double ** dmalloc_2d(int m, int n) {
    if (m <= 0 || n <= 0) return NULL;
    double **A = malloc(m * sizeof(double *));
    if (A == NULL) return NULL;
    A[0] = malloc(m*n*sizeof(double));
    if (A[0] == NULL) {
        free(A);
        return NULL;
 }
    for (int i = 1; i < m; i++)
        A[i] = A[0] + i * n;
    return A;}
// Free it again
void dfree_2d(double **A) {
 free(A[0]);
 free(A);
};


void matmult_nat(int m, int n, int k, double **A, double **B, double **C){
for (int i=0;i<m;i++) {
    for (int p=0;p<n;p++) { 
        for (int j=0;j<k;j++){
            C[i][p]=A[i][j]*B[j][p];
        }
    }
}
}

void matmult_mkn(int m, int n, int k, double **A, double **B, double **C){
for (int i=0;i<m;i++) {
    for (int j=0;j<k;j++){
        for (int p=0;p<n;p++) { 
            C[i][p]=A[i][j]*B[j][p];
        }
    }
}
}

void matmult_nkm(int m, int n, int k, double **A, double **B, double **C){
for (int p=0;p<n;p++) {
    for (int j=0;j<k;j++){
        for (int i=0;i<m;i++) { 
            C[i][p]=A[i][j]*B[j][p];
        }
    }
}
}

void matmult_nmk(int m, int n, int k, double **A, double **B, double **C){
for (int p=0;p<n;p++) {
    for (int i=0;i<m;i++){
        for (int j=0;j<k;j++){
            C[i][p]=A[i][j]*B[j][p];
        }
    }
}
}

void matmult_kmn(int m, int n, int k, double **A, double **B, double **C){
for (int j=0;j<k;j++) {
    for (int i=0;i<m;i++){
        for (int p=0;p<n;p++){
            C[i][p]=A[i][j]*B[j][p];
        }
    }
}
}

void matmult_knm(int m, int n, int k, double **A, double **B, double **C){
for (int j=0;j<k;j++) {
    for (int p=0;p<n;p++){
        for (int i=0;i<m;i++){
            C[i][p]=A[i][j]*B[j][p];
        }
    }
}
}


void main()
{

int k=2;
int n=2;
int m=2;
int i;
int p;
int j;



// Allocating
double **m1=dmalloc_2d(m,k);
double **m2=dmalloc_2d(k,n);
double **m3=dmalloc_2d(m,n);

// Fill with numbers
for (i=0;i<m;i++){
    for (p=0;p<k;p++){
        m1[i][p]=1;
    }
}
for (i=0;i<k;i++){
    for (p=0;p<n;p++){
        m2[i][p]=12000;
    }
}
m2[1][1]=100;
for (i=0;i<m;i++){
    for (p=0;p<n;p++){
        m3[i][p]=0;
    }
}

for (i=0;i<m;i++){
    for (p=0;p<n;p++){
        printf("%f",m3[i][p]);
    }
}


// Deallocating
dfree_2d(m1);
dfree_2d(m2);
dfree_2d(m3);
}