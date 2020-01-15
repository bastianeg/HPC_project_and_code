/* gauss_seidel.c - Poisson problem in 3d
 *
 */
#include <math.h>

void
gauss_seidel(int tol,double ***U, double ***F, int N, int iter_max) {
    double deltasq=((2/N)^2);
    double d = tol+10; //inf
    int iter = 0;
    int tempval;
    double squarenorm;
    while((tol<d) && (iter_max >= iter))
    {
        squarenorm=0;
        // from  i to j to k
        // for i
        for (int i = 1; i<(N-1); i++){
            //for j
            for (int j = 1; j<(N-1); j++){
                //for k
                for (int k = 1; k<(N-1); k++){

                    // U = 1/6 * (sum of us)
                    tempval = U[i][j][k] ;
                    U[i][j][k] = (1/6)*(U[i-1][j][k]+U[i+1][j][k]+U[i][j-1][k]+U[i][j+1][k]+U[i][j][k-1]+U[i][j][k+1]+deltasq*F[i][j][k]);
                    squarenorm += (U[i][j][k]-tempval)*(U[i][j][k]-tempval);
                }
            }
        }
    // norm calc
    d = sqrt(squarenorm);

    // update iteration and Uold
    iter ++;
}
}

