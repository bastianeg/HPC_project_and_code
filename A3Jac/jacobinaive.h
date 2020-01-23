/* jacobi.h - Poisson problem 
 *
 * $Id: jacobi.h,v 1.1 2006/09/28 10:12:58 bd Exp bd $
 */

#ifndef __JACOBINAIVE
#define __JACOBINAIVE

void jacobinaive(double **U, double *F, double **Uold, int N, int iter_max,double **tmp);

#endif
