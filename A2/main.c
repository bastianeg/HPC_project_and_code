/* main.c - Poisson problem in 3D
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "alloc3d.h"
#include "print.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _JACOBI
#include "jacobi.h"
void jacobi(double ***U, double ***F, double ***Uold, int N, int iter_max, double tol);
#endif

#ifdef _GAUSS_SEIDEL
#include "gauss_seidel.h"
#endif

#ifdef _JACOBI_PAR
#include "jacobi_par.h"
#endif

#ifdef _GAUSS_SEIDEL_PAR
#include "gauss_seidel_par.h"
#endif

#define N_DEFAULT 100


void init_data(int N, double ***U, double ***F, double start_T){

    // Initialize U leveraging first touch
    double x, y, z;
    #pragma omp parallel for schedule(static)
    for (int i = 0; i<=(N+1); i++){
        for(int j = 0; j<=(N+1); j++){
            for (int k = 0; k<=(N+1); k++){
                U[i][j][k] = start_T;
            }
        }
    }

    //fill in boundaries of U
    for (int i = 0; i<=(N+1); i++){
        for(int j = 0; j<=(N+1); j++){
                U[i][N+1][j] = 20;
                U[i][0][j] = 0;
                U[N+1][i][j] = 20;
                U[0][i][j] = 20;
                U[i][j][N+1] = 20;
                U[i][j][0] = 20;
        }
    }

    // make F
    #pragma omp parallel for schedule(static)
    for (int i = 0; i<=N+1; i++){
        x = ((2*i)/(double) (N+1))-1;
        for(int j = 0; j<N+1; j++){
            y = (2*j)/(double) (N+1)-1;
            for (int k = 0; k<N+1; k++){
                z = (2*k)/(double) (N+1)-1;
                //then check conditions
                if ((-1.0 <= x) && (x <= -(3/8.0)) && (-1 <= y) && (y <= (-1/2.0)) && ((-2/3.0) <= z) && (z <= 0)){
                    F[i][j][k] = 200;
                }
                else {
                    F[i][j][k] = 0;
                }
            }
        }
    }

}


int
main(int argc, char *argv[]) {

    int 	N = N_DEFAULT;
    int 	iter_max = 1000;
    double	tolerance;
    double	start_T;
    int		output_type = 0;
    char	*output_prefix = "poisson_res";
    char        *output_ext    = "";
    char	output_filename[FILENAME_MAX];
    double 	***u = NULL;
    double 	***u_old = NULL;
    double  ***f = NULL;


    /* get the paramters from the command line */
    N         = atoi(argv[1]);	// grid size
    iter_max  = atoi(argv[2]);  // max. no. of iterations
    tolerance = atof(argv[3]);  // tolerance
    start_T   = atof(argv[4]);  // start T for all inner grid points
    if (argc == 6) {
	output_type = atoi(argv[5]);  // ouput type
    }

    // allocate memory
    if ( (u = d_malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array u: allocation failed");
        exit(-1);
    }
    if ( (u_old = d_malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array u: allocation failed");
        exit(-1);
    }
    if ( (f = d_malloc_3d(N+2, N+2, N+2)) == NULL ) {
        perror("array u: allocation failed");
        exit(-1);
    }

    /////////////////
    init_data(N, u, f, start_T);

    //--->> iterations
    #ifdef _JACOBI
    jacobi(u, f, u_old, N, iter_max, tolerance);
    #endif

    #ifdef _GAUSS_SEIDEL
    gauss_seidel(u, f, N, iter_max, tolerance);
    #endif

    #ifdef _JACOBI_PAR
    jacobi_par(u, f, u_old, N, iter_max, tolerance);
    #endif

    #ifdef _GAUSS_SEIDEL_PAR
    gauss_seidel_par(u, f, N, iter_max);
    #endif

    // dump  results if wanted
    switch(output_type) {
	case 0:
	    // no output at all
	    break;
	case 3:
	    output_ext = ".bin";
	    sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write binary dump to %s: ", output_filename);
	    print_binary(output_filename, N+2, u);
	    break;
	case 4:
	    output_ext = ".vtk";
	    sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write VTK file to %s: ", output_filename);
	    print_vtk(output_filename, N+2, u);
	    break;
	default:
	    fprintf(stderr, "Non-supported output type!\n");
	    break;
    }

    // de-allocate memory
    free(u);
    free(u_old);
    free(f);

    return(0);
}
