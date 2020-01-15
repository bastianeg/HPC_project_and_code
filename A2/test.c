 
#include <stdio.h>
#include <stdlib.h>
#include "alloc3d.h"

#define Uinit 3


double ***
d_malloc_3d(int m, int n, int k) {

    if (m <= 0 || n <= 0 || k <= 0)
        return NULL;

    double ***array3D = malloc(m * sizeof(double **) +
                               m * n * sizeof(double *) +
                               m * n * k * sizeof(double));
    if (array3D == NULL) {
        return NULL;
    }

    for(int i = 0; i < m; i++) {
        array3D[i] = (double **) array3D + m + i * n ;
    }

    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            array3D[i][j] = (double *) array3D + m + m * n + i * n * k + j * k;
        }
    }

    return array3D;
}





int
main(int argc, char *argv[]) {
    
    int N = 4;
    double ***U;
    double ***F;
    double x, y, z;


    U = alloc_3d(N+2, N+2, N+2);
    F = alloc_3d(N+2, N+2, N+2);





    // Fill in U
    for (int i = 1; i<=N; i++){
        for(int j = 1; j<N; j++){
            for (int k = 1; k<N; k++){
                U[i][j][k] = Uinit; //could change this later
            }
        }
    }

    //fill in boundaries of U
    for (int i = 0; i<=N+1; i++){
        for(int j = 0; j<N+1; j++){
                U[i][1][j] = 20;
                U[i][-1][j] = 0;
                U[1][i][j] = 20;
                U[-1][i][j] = 20;
                U[i][j][1] = 20;
                U[i][j][-1] = 20;
        }
    }

    // make F
    for (int i = 0; i<=N+1; i++){

        x = (2*i)/(N+1)-1;

        for(int j = 0; j<N+1; j++){

            y = (2*j)/(N+1)-1;

            for (int k = 0; k<N+1; k++){

                z = (2*k)/(N+1)-1;


                //then check conditions 
                if ((-1 <= x <= (3/8)) && (-1 <= y <= (-1/2)) && ((-2/3) <= z <= 0)){
                    U[i][j][k] = 200;
                }
                else {
                    U[i][j][k] = 0;
                }
                
            }
        }
    }


}
