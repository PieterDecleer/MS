#ifndef _SUBGRID_H              
#define _SUBGRID_H            

#include "fdtd3d.h"
#include <math.h>


#define tol 1e-9     // grid tolerance



// declare static LAPACK functions (use function from external Lapack library which should be linked by gcc compiler) 
static long dgtsv(long N, long NRHS, double *DL, double *D, double *DU, double *B,long LDB){
  extern void dgtsv_(const long *Np, const long *NRHSp, double *DL, double *D, double *DU, double *B, const long *LDBp, long *INFOp);
  long info;
  dgtsv_(&N, &NRHS, DL, D, DU, B, &LDB, &info);
  return info;
}

// coupling interpolation and restriction weights
static double *Wix, *Wiy, *Wiz, *wrx, *Wry, *Wrz;   



// Macros for accessing static variables of Subgrid struct s
#define Wix(I,J,K)       (*g).sg[s].wix[I,J,K]
#define Wiy(I,J,K)       (*g).sg[s].wiy[I,J,K]
#define Wiz(I,J,K)       (*g).sg[s].wix[I,J,K]
#define Wrx(I,J,K)       (*g).sg[s].wrx[I,J,K]
#define Wry(I,J,K)       (*g).sg[s].wry[I,J,K]
#define Wrz(I,J,K)       (*g).sg[s].wrx[I,J,K]


#endif
