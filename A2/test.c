 
#include <stdio.h>
#include <stdlib.h>
#include "alloc3d.h"

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
main(void) {
    
    int N = 4;
    double Uinit=3; //could change this later
    double ***U;
    double ***F;
    double x, y, z;

    U = d_malloc_3d(N+2, N+2, N+2);
    F = d_malloc_3d(N+2, N+2, N+2);

    // Fill in U
    for (int i = 1; i<=N; i++){
        for(int j = 1; j<=N; j++){
            for (int k = 1; k<=N; k++){
                U[i][j][k] = Uinit; 
            }
        }
    }

    //fill in boundaries of U
    
    for (int i = 0; i<=N+1; i++){
        for(int j = 0; j<=N+1; j++){
                U[i][N+1][j] = 20.0;
                U[i][0][j] = 0.0;
                U[N+1][i][j] = 20.0;
                U[0][i][j] = 20.0;
                U[i][j][N+1] = 20.0;
                U[i][j][0] = 20.0;
        }
    }

    // make F
    for (int i = 0; i<=N+1; i++){
        x = ((2*i)/(double) (N+1))-1;
        for(int j = 0; j<N+1; j++){
            y = (2*j)/(double) (N+1)-1;
            for (int k = 0; k<N+1; k++){
                z = (2*k)/(double) (N+1)-1;
                //then check conditions 
                if ((-1.0 <= x) && (x <= -(3/8.0)) && (-1 <= y) && (y <= (-1/2.0)) && ((-2/3.0) <= z) && (z <= 0)){
                    F[i][j][k] = 200.0;
                }
                else {
                    F[i][j][k] = 0.0;
                }
            }
        }
    }


/*
//print F

    // make F
    int k = 0;
    for (int i = 0; i<=N+1; i++){
        for(int j = 0; j<=N+1; j++){
                printf("  %lf  ", F[i][j][k]);
        }
        printf("\n");   
    }
    */
}
