/* gauss_seidel.c - Poisson problem in 3d
 *
 */
#include <math.h>

void
gauss_seidel(tol,double ***U) {
    while((tol<d) && (iter_max >= iter))
    {
        // from  i to j to k
        // for i
        for (int i = 1; i<(N-1); i++){
            //for j
            for (int j = 1; j<(N-1; j++){
                //for k
                for (int k = 2; k<(N-1; k++){

                    // U = 1/6 * (sum of us)
                    U[i][j][k] = (1/6)*(U[i-1][j][k]+U[i+1][j][k]+U[i][j-1][k]+U[i][j+1][k]+U[i][j][k-1]+U[i][j][k+1]+F[i][j][k]);

                    squarenorm += (U[i][j][k]-Uold[i][j][k])*(U[i][j][k]-Uold[i][j][k]);
                }
            }
        }
}

