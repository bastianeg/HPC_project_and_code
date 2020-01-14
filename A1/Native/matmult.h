#ifndef __MATMULT_LIB_H
#define __MATMULT_LIB_H

void matmult_nat(int m, int n, int k, double **A, double **B, double **C){
for (int i=0;i<m;i++){
    for (int p=0;p<n;p++){
        C[i][p]=0;
    }
}
for (int i=0;i<m;i++) {
    for (int p=0;p<n;p++) { 
        for (int j=0;j<k;j++){
            C[i][p]+=A[i][j]*B[j][p];
        }
    }
}
}

void matmult_mkn(int m, int n, int k, double **A, double **B, double **C){
for (int i=0;i<m;i++){
    for (int p=0;p<n;p++){
        C[i][p]=0;
    }
}
for (int i=0;i<m;i++) {
    for (int j=0;j<k;j++){
        for (int p=0;p<n;p++) { 
            C[i][p]+=A[i][j]*B[j][p];
        }
    }
}
}

void matmult_nkm(int m, int n, int k, double **A, double **B, double **C){
for (int i=0;i<m;i++){
    for (int p=0;p<n;p++){
        C[i][p]=0;
    }
}
for (int p=0;p<n;p++) {
    for (int j=0;j<k;j++){
        for (int i=0;i<m;i++) { 
            C[i][p]+=A[i][j]*B[j][p];
        }
    }
}
}

void matmult_nmk(int m, int n, int k, double **A, double **B, double **C){
for (int i=0;i<m;i++){
    for (int p=0;p<n;p++){
        C[i][p]=0;
    }
}
for (int p=0;p<n;p++) {
    for (int i=0;i<m;i++){
        for (int j=0;j<k;j++){
            C[i][p]+=A[i][j]*B[j][p];
        }
    }
}
}

void matmult_kmn(int m, int n, int k, double **A, double **B, double **C){
for (int i=0;i<m;i++){
    for (int p=0;p<n;p++){
        C[i][p]=0;
    }
}
for (int j=0;j<k;j++) {
    for (int i=0;i<m;i++){
        for (int p=0;p<n;p++){
            C[i][p]+=A[i][j]*B[j][p];
        }
    }
}
}

void matmult_knm(int m, int n, int k, double **A, double **B, double **C){
for (int i=0;i<m;i++){
    for (int p=0;p<n;p++){
        C[i][p]=0;
    }
}

for (int j=0;j<k;j++) {
    for (int p=0;p<n;p++){
        for (int i=0;i<m;i++){
            C[i][p]+=A[i][j]*B[j][p];
        }
    }
}
}

#endif