/* jacobi.c - Poisson problem in 3d
 *
 */
#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>
#include <stdio.h>

void
jacobi(double ***U, double ***F, double ***Uold, int N, int iter_max, double tol) {
    double ts, te; // for timing
    double deltasq = 4.0/((double) N * (double) N);
    //define norm and max_iter and Uold and iter and threshold
    int iter = 0;
    double onesixth = 1.0/6.0;
    double d = tol*N*N*N+10; //inf

    // update Uold = U
    for (int i = 0; i<(N+2); i++){
        for (int j = 0; j<(N+2); j++){
            for (int k = 0; k<(N+2); k++){
                Uold[i][j][k] = U[i][j][k];
                F[i][j][k] *= deltasq;
            }
        }
    }
    ts = omp_get_wtime();
    //while condition is not satisfied
    while((d/sqrt(N*N*N)>tol) && (iter < iter_max))
    {
        // start wallclock timer
        
        d = 0.0;

        // from  i to j to k
        // for i
        for (int i = 1; i<(N+1); i++){
            //for j
            for (int j = 1; j<(N+1); j++){
                //for k
                for (int k = 1; k<(N+1); k++){

                    // U = 1/6 * (sum of us +Delta^2 f)
                    U[i][j][k] = onesixth*(Uold[i-1][j][k]+Uold[i+1][j][k]+Uold[i][j-1][k]+\
                    Uold[i][j+1][k]+Uold[i][j][k-1]+Uold[i][j][k+1]+F[i][j][k]);
                    
                    // frobenius norm
                    d += (U[i][j][k]-Uold[i][j][k])*(U[i][j][k]-Uold[i][j][k]);
                }
            }
        }

        // norm calc
        d = sqrt(d);

        // update iteration and Uold
        iter ++;

        for (int i = 0; i<(N+2); i++){
            for (int j = 0; j<(N+2); j++){
                for (int k = 0; k<(N+2); k++){
                    Uold[i][j][k] = U[i][j][k];
                }
            }
        }
    }
    te = omp_get_wtime() - ts;
    
    printf("Number of iterations: %d\n", iter);
    printf("Norm: %lf\n", d);
    printf("Elapsed time: %lf\n", te);
    printf("Iterations per second: %lf\n", iter/te);
    //printf("%.5lf, %.5lf\n", te, 1e-6*11*N*N*N*iter/te);
}
