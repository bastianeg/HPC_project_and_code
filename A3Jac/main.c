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

#ifdef _JACOBINAIVE
#include "jacobinaive.h"
#endif

#ifdef _JACOBISEQ
#include "jacobiseq.h"
#endif

#ifdef _JACOBIMULTI
#include "jacobi_multi.h"
#endif

#ifdef _JACOBITOL
#include "jacobi_tol.h"
#endif

#define N_DEFAULT 1


void init_data(int N, double *U, double *F, double start_T){
    // Initialize U leveraging first touch
    double x, y, z;
    for (int i = 0; i<=(N+1); i++){
        for(int j = 0; j<=(N+1); j++){
            for (int k = 0; k<=(N+1); k++){
                U[i+j+k] = start_T;
            }
        }
    }

    //fill in boundaries of U
    for (int i = 0; i<=(N+1); i++){
        for(int j = 0; j<=(N+1); j++){
                U[i+(N+2)*(N+1)+(N+2)*(N+2)*j] = 20;
                U[i+(N+2)*(N+2)*j] = 0;
                U[N+1+(N+2)*i+(N+2)*(N+2)*j] = 20;
                U[(N+2)*i+(N+2)*(N+2)*j] = 20;
                U[i+j*(N+2)+(N+2)*(N+2)*(N+1)] = 20;
                U[i+j*(N+2)] = 20;
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
                    F[i+(N+2)*j+(N+2)*(N+2)*k] = 200;
                }
                else {
                    F[i+(N+2)*j+(N+2)*(N+2)*k] = 0;
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
    double 	*u = NULL;
    double 	*u_old = NULL;
    double  *f = NULL;


    /* get the paramters from the command line */
    N         = atoi(argv[1]);	// grid size
    iter_max  = atoi(argv[2]);  // max. no. of iterations
    tolerance = atof(argv[3]);  // tolerance
    start_T   = atof(argv[4]);  // start T for all inner grid points
    if (argc == 6) {
	output_type = atoi(argv[5]);  // ouput type
    }

    // allocate memory (properly!)
    if ( u = malloc((N+2)*(N+2)*(N+2)*sizeof(double)) == NULL ) {
        perror("array u: allocation failed");
        exit(-1);
    }
    if ( u_old = malloc((N+2)*(N+2)*(N+2)*sizeof(double)) == NULL ) {
        perror("array u: allocation failed");
        exit(-1);
    }
    if ( f = malloc((N+2)*(N+2)*(N+2)*sizeof(double)) == NULL ) {
        perror("array u: allocation failed");
        exit(-1);
    }

    /////////////////
    init_data(N, u, f, start_T);

    //allocate memory on GPU
    double* D_u;
    double* D_u_old;
    double* D_f;
    cudaMalloc((void**) &D_u, (N+2)*(N+2)*(N+2)*sizeof(double));
    cudaMalloc((void**) &D_u_old, (N+2)*(N+2)*(N+2)*sizeof(double));
    cudaMalloc((void**) &D_f, (N+2)*(N+2)*(N+2)*sizeof(double));

    //move u to GPU
    cudaMemcpy(D_u, u, (N+2)*(N+2)*(N+2)*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(D_f, f, (N+2)*(N+2)*(N+2)*sizeof(double), cudaMemcpyHostToDevice);

    //--->> iterations
    #ifdef _JACOBISEQ
    jacobiseq(u, f, u_old, N, iter_max);
    #endif

    #ifdef _GAUSSNAIVE
    gaussnaive(u, f, N, iter_max);
    #endif
    
    #ifdef _JACOBIMULT
    cudaSetDevice(0);
    double *d0_U;
    cudaMalloc((void**)&d0_U, (N+2)*(N+2)*(N+2)/2*sizeof(double));
    cudaMemcpy(d0_U, u, (N+2)*(N+2)*(N+2)/2*sizeof(double), cudaMemcpyHostToDevice);

    cudaSetDevice(1);
    double *d1_U;
    cudaMalloc((void**)&d1_U, A_size/2);
    cudaMemcpy(d1_U, u + (N+2)*(N+2)*(N+2)/2, (N+2)*(N+2)*(N+2)/2*sizeof(double), cudaMemcpyHostToDevice);

    jacobimulti(d0_U, d1_U, f, u_old, N, iter_max);
    #endif

    #ifdef _JACOBITOL
    jacobitol(u, f, N, iter_max,tolerance);
    #endif

    //move u back to host
    cudaMemcpy(u, D_u, (N+2)*(N+2)*(N+2)*sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(D_u);
    cudaFree(D_u_old);
    cudaFree(D_f);

    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            printf("%.2lf ",u[i*N+j]);
        }
        printf("\n");
    }
    printf("on the CPU\n");

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
	// case 4:
	//     output_ext = ".vtk";
	//     sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	//     fprintf(stderr, "Write VTK file to %s: ", output_filename);
	//     print_vtk(output_filename, N+2, u);
	//     break;
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
