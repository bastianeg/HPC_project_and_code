#include <stdio.h>
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

void matmult_blk(int m,int n,int k,double **A,double **B,double **C, int bs){

    int i,j,l, iend, jend, lend;
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
    //From our experiments, mkn is the fastest across all problem sizes.
    for(bi=0; bi<mb; bi++){
        for(bl=0; bl<kb; bl++){
            for(bj=0; bj<nb; bj++){

                //looping within a block
                //we make sure the i,l,j is the actual i,l,j
                //to account non multiples of bs we use MIN
                iend = MIN((bi+1)*bs,m);
                for(i=bi*bs; i<iend; i++){
                    lend = MIN((bl+1)*bs,k);
                    for(l=bl*bs; l<lend; l++){
                        jend = MIN((bj+1)*bs,n);
                        for(j=bj*bs; j<jend; j++){
                            C[i][j] += A[i][l]*B[l][j];
                        }
                    }
                }
            }
        }
    }

}
