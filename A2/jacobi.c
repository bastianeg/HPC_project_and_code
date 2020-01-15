/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <math.h>



void
jacobi(double ***U, double ***F, double ***Uold, int N, int iter_max, double *tol) {
    // fill in your code here

    //define norm and max_iter and Uold and iter and threshold
    double U1, U2, U3, U4, U5, U6, b, squarenorm;
    int iter = 0;
    double d = tol+10; //inf
    
    Uold = U;
    
    //while condition is not satisfied
    while((tol<d) || (iter_max >= iter))
    {
        // from  i to j to k
        // for i
        for (int i = 1; i<(N-1); i++){
            //for j
            for (int j = 1; j<(N-1; j++){
                //for k
                for (int k = 2; k<(N-1; k++){

                    //update all Us
                    U1 = Uold[i-1][j][k];
                    U2 = Uold[i+1][j][k];
                    U3 = Uold[i][j-1][k];
                    U4 = Uold[i][j+1][k];
                    U5 = Uold[i][j][k-1];
                    U6 = Uold[i][j][k+1];
                    b = F[i][j][k];

                    // U = 1/6 * (sum of us)
                    U[i][j][k] = (1/6)*(U1+U2+U3+U4+U5+U6+b);

                    squarenorm += (U[i][j][k]-Uold[i][j][k])*(U[i][j][k]-Uold[i][j][k]);
                }
            }
        }
    
    // norm calc
    d = sqrt(squarenorm);

    // update iteration and Uold
    iter ++;
    Uold = U;
    }      
}
