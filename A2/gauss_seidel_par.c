/* gauss_seidel.c - Poisson problem in 3d
 *
 */
#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>
// #include <stdio.h>

void
gauss_seidel_par(double ***U, double ***F, int N, int iter_max) {
    double deltasq = 4.0/((double) N * (double) N);
    double onesixth = 1.0/6.0;
    int iter = 0;
    int tempval;
    double squarenorm;
    while(iter < iter_max)
    {
        #pragma omp parallel default(none) shared(N,onesixth,U,F,deltasq)
        {
            
            #pragma omp for ordered(1)
            for (int i = 1; i<(N+1); i++){
                for (int j = 1; j<(N+1); j++){
                    for (int k = 1; k<(N+1); k++){
                        
                        #pragma omp ordered depend(sink:i-1,j-1,k-1)
                        U[i][j][k] = onesixth*(U[i-1][j][k]+U[i+1][j][k]+U[i][j-1][k]+U[i][j+1][k]+U[i][j][k-1]+U[i][j][k+1]+deltasq*F[i][j][k]);
                        #pragma omp ordered depend(source)
                    }
                }
            }
        }//ending parralel region
    // norm calc

    // update iteration and Uold
    iter ++;
}
}
