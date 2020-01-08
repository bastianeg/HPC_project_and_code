#include <stdio.h>
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

void matmult_blk(int m,int n,int k,double **A,double **B,double **C, int bs){

    int i,j,l;
    int bi,bj,bl;

    // The number of blocks is ceil of N/bs
    int mb = m/bs + (int)(m%bs!=0);
    int nb = n/bs + (int)(n%bs!=0);
    int kb = k/bs + (int)(k%bs!=0);

    //initialize C as zeros
    for(i=0; i<m; i++)
        for(j=0; j<n; j++)
            C[i][j] = 0.0;

    //looping through the blocks
    for(bi=0; bi<mb; bi++){
        for(bj=0; bj<nb; bj++){
            for(bl=0; bl<kb; bl++){
                
                //looping within a block
                //we make sure the i,j,l is the actual i,j,l
                //to account non multiples of bs we use MIN
                for(i=bi*bs; i<MIN((bi+1)*bs,m); i++){
                    for(j=bj*bs; j<MIN((bj+1)*bs,n); j++){
                        for(l=bl*bs; l<MIN((bl+1)*bs,k); l++){
                            C[i][j] += A[i][l]*B[l][j];
                        }
                    }
                }
            }
        }
    }

}