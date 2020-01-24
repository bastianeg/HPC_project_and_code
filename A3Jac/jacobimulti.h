/* jacobi.h - Poisson problem 
 *
 * $Id: jacobi.h,v 1.1 2006/09/28 10:12:58 bd Exp bd $
 */

#ifndef __JACOBIMULTI
#define __JACOBIMULTI

void jacobimulti(double* D0U,double* D1U, double* D0F, double* D1F, double* D0Uold, double* D1old, int N, int iter_max);

#endif
