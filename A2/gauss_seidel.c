/* gauss_seidel.c - Poisson problem in 3d
 *
 */
#include <math.h>

void
gauss_seidel(double ***U, double ***F, int N, int iter_max,double tol) {
    double deltasq = 4.0/((double) N * (double) N);
    double onesixth = 1.0/6.0;
    double d = tol+10.0; //inf
    int iter = 0; 
    int tempval;
    double squarenorm;
    while((tol<d) || (iter_max >= iter))
    {
        d=0;
        // from  i to j to k
        // for i
        for (int i = 1; i<(N+1); i++){
            //for j
            for (int j = 1; j<(N+1); j++){
                //for k
                for (int k = 1; k<(N+1); k++){

                    // U = 1/6 * (sum of us)
                    tempval = U[i][j][k];
                    U[i][j][k] = onesixth*(U[i-1][j][k]+U[i+1][j][k]+U[i][j-1][k]+U[i][j+1][k]+U[i][j][k-1]+U[i][j][k+1]+deltasq*F[i][j][k]);
                    d += (U[i][j][k]-tempval)*(U[i][j][k]-tempval);
                }
            }
        }
    // norm calc
    d = sqrt(d);

    // update iteration and Uold
    iter ++;
    printf("completed iteration %d\n",iter);
}
}

