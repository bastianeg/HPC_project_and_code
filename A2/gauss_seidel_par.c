/* gauss_seidel.c - Poisson problem in 3d
 *
 */
#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>
#include <stdio.h>

void
gauss_seidel_par(double ***U, double ***F, int N, int iter_max) {
    // parameters and predifinitions
    double deltasq = 4.0/((double) N * (double) N);
    double onesixth = 1.0/6.0;
    int tempval;
    double squarenorm;
    double ts = omp_get_wtime();
    // starting parralel region. Each iteration has a barrier
    #pragma omp parallel default(none) shared(N,onesixth,U,F,deltasq,iter_max)
    {
        // main iteration loop
        for(int iter = 0; iter < iter_max; iter++)
        {
            //starting doacross loop for the 3-nested for loop.
            #pragma omp for ordered(2)
            for (int i = 1; i<(N+1); i++){
                for (int j = 1; j<(N+1); j++){
                    //depending on the two outer dimensions such that we are doing one line in each process.
                    #pragma omp ordered depend(sink:i-1,j) depend(sink:i,j-1)
                    for (int k = 1; k<(N+1); k++){             
                        U[i][j][k] = onesixth*(U[i-1][j][k]+U[i+1][j][k]+U[i][j-1][k]+U[i][j+1][k]+U[i][j][k-1]+U[i][j][k+1]+deltasq*F[i][j][k]);
                    }
                    #pragma omp ordered depend(source)
                }
            }
            
        }
    }//ending parralel region
    double te = omp_get_wtime() - ts;
    printf("%.5lf, %.5lf\n", te, 1e-6*8*N*N*N*iter_max/te);
    printf("Number of iterations: %d\n", iter_max);
    printf("Elapsed time: %lf\n", te);
    printf("Iterations per second: %lf\n", iter_max/te);
}
