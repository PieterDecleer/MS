#ifndef _ADHIE_H              
#define _ADHIE_H            

#include "fdtd3d.h"
                                             

struct Slice{ 
        int i, ne, nh;     // specifies the range of cells that is implicitized: for example, Dz2(i) to Dz2(i+ne-1) is implicitized for the update of Ey            
	double *qe, *qh;   // ADHIE update coefficients for the tridiagonal solver                     
};   
typedef struct Slice Slice;  
                                          

// macros for accessing auxiliary unknowns of the Slice structs                                    
#define Ix(S)           slx[S].i 
#define Nye(S)          slx[S].ne
#define Nyh(S)          slx[S].nh 
#define Qex(S,I,J,K,L)  slx[S].qe[(((I)*Nye(S)+J)*(Nz-1)+K)*3+L] 
#define Qhx(S,I,J,K,L)  slx[S].qh[(((I)*Nyh(S)+J)*Nz+K)*3+L] 

#define Iy(S)           sly[S].i 
#define Nze(S)          sly[S].ne 
#define Nzh(S)          sly[S].nh 
#define Qey(S,I,J,K,L)  sly[S].qe[(((I)*Ny+J)*Nze(S)+K)*3+L]   
#define Qhy(S,I,J,K,L)  sly[S].qh[(((I)*(Ny-1)+J)*Nzh(S)+K)*3+L] 

#define Iz(S)           slz[S].i 
#define Nxe(S)          slz[S].ne 
#define Nxh(S)          slz[S].nh 
#define Qez(S,I,J,K,L)  slz[S].qe[(((I)*(Ny-1)+J)*Nz+K)*3+L]                            
#define Qhz(S,I,J,K,L)  slz[S].qh[(((I)*Ny+J)*(Nz-1)+K)*3+L] 


// declare static variables
static int Sx, Sy, Sz;               // the number of slices or, equivalently, the number tridiagonal blocks 
static Slice *slx, *sly, *slz; 



#endif  





