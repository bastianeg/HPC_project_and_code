 
#include <stdio.h>
#include <stdlib.h>
#include "alloc3d.h"

#define Uinit 3



int
main(int argc, char *argv[]) {
    
    int N = 4;
    double ***U;
    double ***F;
    double x, y, z;


    U = d_malloc_3d(N+2, N+2, N+2);
    F = d_malloc_3d(N+2, N+2, N+2);





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



//print F

    // make F
    int k = 0;
    for (int i = 0; i<=N+1; i++){
        for(int j = 0; j<N+1; j++){
                    printf("%ld ", U[i][j][k]);
        }
        printf("\n");
    }
}
