/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>



void
jacobi(double ***U, double ***F, double ***Uold, int N, int iter_max, double tol) {
    
    double deltasq = 4.0/(double) N*N;
    //define norm and max_iter and Uold and iter and threshold
    double U1, U2, U3, U4, U5, U6;
    int iter = 0;
    double onesixth = 1.0/6.0;
    double d = tol+10; //inf
    for (int i = 0; i<(N+2); i++){
        for (int j = 0; j<(N+2); j++){
            for (int k = 0; k<(N+2); k++){
                Uold[i][j][k] = U[i][j][k];
            }
        }
    } 
    
    //while condition is not satisfied
    while((tol<d) || (iter_max >= iter))
    {
        d = 0.0;

        // from  i to j to k
        // for i
        for (int i = 1; i<(N-1); i++){
            //for j
            for (int j = 1; j<(N-1); j++){
                //for k
                for (int k = 1; k<(N-1); k++){

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
}
