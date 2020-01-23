/* jacobi.h - Poisson problem 
 *
 * $Id: jacobi.h,v 1.1 2006/09/28 10:12:58 bd Exp bd $
 */

#ifndef _JACOBISEQ
#define _JACOBISEQ

void jacobiseq(double ***U, double ***F, double ***Uold, int N, int iter_max, double tol);

__global__ void updmat(int N, double* U, double* Uold, int i, int j, int k);

__global__ void jacgpu(int N, double* A, double* b, double* onesixth, int i, int j, int k);

__global__ void initmat(int N, double* U, double* Uold, double* F,double deltasq, int i, int j, int k);


#endif
