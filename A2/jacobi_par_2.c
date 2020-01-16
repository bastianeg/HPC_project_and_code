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
    double U1, U2, U3, U4, U5, U6;
    double onesixth = 1.0/6.0;
    double d = tol+10; //initialize norm to inf

    ts = omp_get_wtime(); // start wallclock timer

    // update Uold = U
    #pragma omp parallel for
    for (int i = 0; i<(N+2); i++){
        for (int j = 0; j<(N+2); j++){
            for (int k = 0; k<(N+2); k++){
                Uold[i][j][k] = U[i][j][k];
            }
        }
    }


     //while condition is not satisfied
     for (int iter = 0; iter < iter_max; iter++){
         if (d < tol) break;
    
        d = 0.0;

        // from  i to j to k
        #pragma omp parallel for reduction(+ : d) default(none) shared(N, U1, U2, U3, U4, U5, U6, Uold, onesixth, deltasq, F, U)
        {
            for (int i = 1; i<(N+1); i++){
                //for j
                for (int j = 1; j<(N+1); j++){
                    //for k
                    for (int k = 1; k<(N+1); k++){

                        //update all Us
                        U1 = Uold[i-1][j][k];
                        U2 = Uold[i+1][j][k];
                        U3 = Uold[i][j-1][k];
                        U4 = Uold[i][j+1][k];
                        U5 = Uold[i][j][k-1];
                        U6 = Uold[i][j][k+1];

                        // U = 1/6 * (sum of us +Delta^2 f)
                        U[i][j][k] = onesixth*(U1+U2+U3+U4+U5+U6+deltasq*F[i][j][k]);

                        d += (U[i][j][k]-Uold[i][j][k])*(U[i][j][k]-Uold[i][j][k]);
                    }
                } 
            }
        } // end parallel

        // norm calc
        d = sqrt(d);

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
    printf("Elapsed time: %.5lf\n", te);
}
