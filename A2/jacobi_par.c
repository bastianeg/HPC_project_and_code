/* jacobi.c - Poisson problem in 3d
 *
 */
#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>
#include <stdio.h>

void
jacobi_par(double ***U, double ***F, double ***Uold, int N, int iter_max, double tol) {
    double ts, te; // for timing
    double deltasq = 4.0/((double) N * (double) N);
    double onesixth = 1.0/6.0;
    double d = tol+10; //initialize norm to inf
    int iter = 0;

    ts = omp_get_wtime(); // start wallclock timer

    // update Uold = U
    #pragma omp parallel for
    for (int i = 0; i<(N+2); i++){
        for (int j = 0; j<(N+2); j++){
            for (int k = 0; k<(N+2); k++){
                Uold[i][j][k] = U[i][j][k];
                F[i][j][k] *= deltasq;
            }
        }
    }


     //while condition is not satisfied
    while((d>tol) && (iter < iter_max)){
    
        d = 0.0;

        // from  i to j to k
        #pragma omp parallel for reduction(+ : d)
        {
            for (int i = 1; i<(N+1); i++){
                //for j
                for (int j = 1; j<(N+1); j++){
                    //for k
                    for (int k = 1; k<(N+1); k++){

                        // U = 1/6 * (sum of us +Delta^2 f)
                        // this should be unrolled factor 4 in k
                        U[i][j][k] = onesixth*(Uold[i-1][j][k]+Uold[i+1][j][k]+Uold[i][j-1][k]+Uold[i][j+1][k]+Uold[i][j][k-1]+Uold[i][j][k+1]+F[i][j][k]);

                        d += (U[i][j][k]-Uold[i][j][k])*(U[i][j][k]-Uold[i][j][k]);
                    }
                } 
            }
        } // end parallel

        // norm calc
        d = sqrt(d);
        iter++;

        // update Uold
        #pragma omp parallel for
        for (int i = 0; i<(N+2); i++){
            for (int j = 0; j<(N+2); j++){
                for (int k = 0; k<(N+2); k++){
                    Uold[i][j][k] = U[i][j][k];
                }
            }
        }
    }
    te = omp_get_wtime() - ts;
    printf("%.5lf, %.5lf\n", te, 1e-9*11*N*N*N*iter/te);
}
