/* main.c - Poisson problem in 3D
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "print.h"
#include <omp.h>


#ifdef _JACOBISEQ
#include "jacobiseq.h"
#endif

#ifdef _JACOBINAIVE
#include "jacobinaive.h"
#endif

#ifdef _JACOBIMULTI
#include "jacobimulti.h"
#endif

#ifdef _JACOBITOL
#include "jacobitol.h"
#endif

#define N_DEFAULT 128

void init_data(int N, double *U, double *F, double start_T){
    // Initialize U leveraging first touch
    double x, y, z;
    int i,j,k;
    for (i = 0; i<=(N+1); i++){
        for(j = 0; j<=(N+1); j++){
            for (k = 0; k<=(N+1); k++){
                U[i*(N+2)*(N+2)+j*(N+2)+k] = start_T;
            }
        }
    }

    //fill in boundaries of U
    for (i = 0; i<=(N+1); i++){
        for(j = 0; j<=(N+1); j++){
                U[i+(N+2)*(N+1)+(N+2)*(N+2)*j] = 20;
                U[i+(N+2)*(N+2)*j] = 0;
                U[N+1+(N+2)*i+(N+2)*(N+2)*j] = 20;
                U[(N+2)*i+(N+2)*(N+2)*j] = 20;
                U[i+j*(N+2)+(N+2)*(N+2)*(N+1)] = 20;
                U[i+j*(N+2)] = 20;
        }
    }

    // make F
    for (i = 0; i<=N+1; i++){
        x = ((2*i)/(double) (N+1))-1;
        for(j = 0; j<N+1; j++){
            y = (2*j)/(double) (N+1)-1;
            for (k = 0; k<N+1; k++){
                z = (2*k)/(double) (N+1)-1;
                //then check conditions
                if ((-1.0 <= x) && (x <= -(3/8.0)) && (-1 <= y) && (y <= (-1/2.0)) && ((-2/3.0) <= z) && (z <= 0)){
                    F[i+(N+2)*j+(N+2)*(N+2)*k] = 200.0;
                }
                else {
                    F[i+(N+2)*j+(N+2)*(N+2)*k] = 0.0;
                }
            }
        }
    }

}


int
main(int argc, char *argv[]) {

    int 	N = N_DEFAULT;
    double	start_T;
    int		output_type = 0;
    double*  u = NULL;
    double*  u_old = NULL;
    double*  f = NULL;
    int i,j,k;
    int iter_max = 100;

    /* get the paramters from the command line */
    N         = atoi(argv[1]);	// grid size
    iter_max  = atoi(argv[2]);  // max. no. of iterations
    start_T   = atof(argv[4]);  // start T for all inner grid points
    if (argc == 6) {
	output_type = atoi(argv[5]);  // ouput type
    }

    // allocate memory
    u = (double*) malloc((N+2)*(N+2)*(N+2)*sizeof(double));
    if ( u == NULL ) {
        perror("array u: allocation failed");
        exit(-1);
    }
    u_old = (double*) malloc((N+2)*(N+2)*(N+2)*sizeof(double));
    if ( u_old == NULL ) {
        perror("array u: allocation failed");
        exit(-1);
    }
    f = (double*) malloc((N+2)*(N+2)*(N+2)*sizeof(double));
    if ( f == NULL ) {
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
    jacobiseq(D_u, D_f, D_u_old, N, iter_max);
    #endif

    #ifdef _JACOBINAIVE
    jacobinaive(D_u, D_f, D_u_old, N, iter_max);
    #endif
    
    #ifdef _JACOBIMULTI
    cudaSetDevice(0);
    double *d0_U;
    double *d0_Uold;
    double *d0_F;
    cudaDeviceEnablePeerAccess(1, 0);
    cudaMalloc((void**)&d0_U, (N+2)*(N+2)*(N+2)/2*sizeof(double));
    cudaMemcpy(d0_U, u, (N+2)*(N+2)*(N+2)/2*sizeof(double), cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d0_Uold, (N+2)*(N+2)*(N+2)/2*sizeof(double));
    cudaMemcpy(d0_Uold, u_old, (N+2)*(N+2)*(N+2)/2*sizeof(double), cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d0_F, (N+2)*(N+2)*(N+2)/2*sizeof(double));
    cudaMemcpy(d0_F, f , (N+2)*(N+2)*(N+2)/2*sizeof(double), cudaMemcpyHostToDevice);


    cudaSetDevice(1);
    double *d1_U;
    double *d1_Uold;
    double *d1_F;
    cudaDeviceEnablePeerAccess(0, 0);
    cudaMalloc((void**)&d1_U, (N+2)*(N+2)*(N+2)/2*sizeof(double));
    cudaMemcpy(d1_U, u + (N+2)*(N+2)*(N+2)/2, (N+2)*(N+2)*(N+2)/2*sizeof(double), cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d1_Uold, (N+2)*(N+2)*(N+2)/2*sizeof(double));
    cudaMemcpy(d1_Uold, u_old + (N+2)*(N+2)*(N+2)/2, (N+2)*(N+2)*(N+2)/2*sizeof(double), cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d1_F, (N+2)*(N+2)*(N+2)/2*sizeof(double));
    cudaMemcpy(d1_F, f + (N+2)*(N+2)*(N+2)/2, (N+2)*(N+2)*(N+2)/2*sizeof(double), cudaMemcpyHostToDevice);

    jacobimulti(d0_U, d1_U, d0_F, d1_F, d0_Uold, d1_Uold, N, iter_max);

    //move back and merge into one array
    cudaMemcpy(u,d0_U,(N+2)*(N+2)*(N+2)/2*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(u+(N+2)*(N+2)*(N+2)/2,d1_U,(N+2)*(N+2)*(N+2)/2*sizeof(double), cudaMemcpyDeviceToHost);

    #endif

    #ifdef _JACOBITOL
    double tolerance = 1.5e-3;
    tolerance = atof(argv[3]);  // tolerance
    double* res;
    cudaMalloc((void**) &res, (N+2)*(N+2)*(N+2)*sizeof(double));
    jacobitol(D_u, D_f, D_u_old, N, iter_max,tolerance,res);
    #endif

    //move u back to host
    #ifndef _JACOBIMULTI
    cudaMemcpy(u, D_u, (N+2)*(N+2)*(N+2)*sizeof(double), cudaMemcpyDeviceToHost);
    #endif

    cudaFree(D_u);
    cudaFree(D_u_old);
    cudaFree(D_f);
    
    for(i = 0; i<2; i++){
        printf("k=%d\n",i-1);
        for(j = 0; j<N+2; j++){
            for(k = 0; k<N+2; k++){
                printf("%6.2lf ",u[i*(N+2)*(N+2)+j*(N+2)+k]);
            }
            printf("\n");
        }
        printf("\n\n");
        
    }
    printf("on the CPU\n");
    

    // dump  results if wanted
    switch(output_type) {
	case 0:
	    // no output at all
	    break;
	// case 3:
	//     output_ext = ".bin";
	//     sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	//     fprintf(stderr, "Write binary dump to %s: ", output_filename);
	//     print_binary(output_filename, N+2, u);
	//     break;
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
