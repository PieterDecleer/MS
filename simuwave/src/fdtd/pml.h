#ifndef _PML_H              
#define _PML_H            

#include "fdtd3d.h"
#include <math.h>


#define m_pg      3.0     
#define Amax(I)   (PyFloat_AS_DOUBLE(PyList_GET_ITEM(py_Amax,I)))  
#define Kmax(I)   (PyFloat_AS_DOUBLE(PyList_GET_ITEM(py_Kmax,I)))  
                                             

struct Pml{ 
	double *e1, *e2, *h1, *h2;         // tangential field components according to the righ-hand rule (e.g. for the x-pml: 1=y and 2=z)  
	double *ae, *ah;                   // update coefficients of the auxiliary equations
	int ntan[2];                       // number of discretization points in the tangential directions
};   
typedef struct Pml Pml;  
                                          

// macros for accessing auxiliary unknowns of the Pml struct                                    
#define Ntan(P,I)     pmls[P].ntan[I]                                  // P is the PML number ({0,1,2,3,4,5}={-x,+x,-y,+y,-z,+z}) 
#define E1(P,I,J,K)   pmls[P].e1[((I)*Ntan(P,0)+J)*(Ntan(P,1)-1)+K]    // I corresponds to the normal direction,  
#define E2(P,I,J,K)   pmls[P].e2[((I)*(Ntan(P,0)-1)+J)*Ntan(P,1)+K]    // J to the first tangential direction,
#define H1(P,I,J,K)   pmls[P].h1[((I)*(Ntan(P,0)-1)+J)*Ntan(P,1)+K]    // and K to the second tangential direction
#define H2(P,I,J,K)   pmls[P].h2[((I)*Ntan(P,0)+J)*(Ntan(P,1)-1)+K]    // according to the right-hand rule
#define Ae(P,I,U)     pmls[P].ae[(I)*2+U]                              
#define Ah(P,I,U)     pmls[P].ah[(I)*2+U]                                

// memory allocation macro
#define my_alloc_pml(P)                                                     \
	my_alloc( pmls[P].e1,  Npml(P)*Ntan(P,0)*(Ntan(P,1)-1) , double );  \
	my_alloc( pmls[P].e2,  Npml(P)*(Ntan(P,0)-1)*Ntan(P,1) , double );  \
	my_alloc( pmls[P].h1,  Npml(P)*(Ntan(P,0)-1)*Ntan(P,1) , double );  \
	my_alloc( pmls[P].h2,  Npml(P)*Ntan(P,0)*(Ntan(P,1)-1) , double );  \
	my_alloc( pmls[P].ae,  Npml(P)*2                       , double );  \
	my_alloc( pmls[P].ah,  Npml(P)*2                       , double );  \


// declare static variables
static Pml pmls[6]; 


#endif  


/* COMMENTS: 
(1) We use the complex-frequency-shifted (CFS) perfectly-matched layer (PML) with recursive convolution (CPML).  
(2) The major E and H fields have the same update equations (with varying update coefficients) throughout the 
     whole 3D space. The auxiliary variables inside the PML are updated in a subroutine and are injected into 
     the major fields (see appendix B of Gedney's book).   
*/
