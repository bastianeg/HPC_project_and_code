/* jacobi.h - Poisson problem 
 *
 * $Id: jacobi.h,v 1.1 2006/09/28 10:12:58 bd Exp bd $
 */

#ifndef __JACOBIMULTI
#define __JACOBIMULTI

void jacobimulti(double *D0U,double* D1U, double *F, double *Uold, int N, int iter_max);

#endif
