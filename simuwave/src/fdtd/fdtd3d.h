#ifndef _FDTD3D_H              
#define _FDTD3D_H            

#include <Python.h>
#include <stdio.h>
#include <stdlib.h>  


// vacuum parameters  
#define Z0 376.730313461
#define c0 299792458.0
#define Hbar	(1.0545718e-34)
#define Qe	(-1.60217662e-19)
#define	Me	(9.10938356e-31)


// memory allocation macro
#define my_alloc(pntr,num,type)                                                           \
	pntr = (type*)calloc(num,sizeof(type));                                           \
	if (!pntr){                                                                       \
		perror("my_alloc");                                                       \
		fprintf(stderr,"Allocation failed for " #pntr ". Terminating... \n");     \
		exit(EXIT_FAILURE);                                                       \
	}   


// Python list macro 
// 'PyFloat_FromDouble' creates a new reference, which is stolen by 'PyList_SetItem' 
// and is garbage-collected once 'Py_DECREF' is applied to 'list'. 
#define my_SetItem(list,ind,item)                                                 \
	int check = PyList_SetItem(list,ind,PyFloat_FromDouble(item));            \
	if (check!=0){ printf("PyList_SetItem failed.\n"); exit(EXIT_FAILURE); }  \



// minimum and maximum macros
# define min(a,b) ((a)<(b) ? (a):(b))
# define max(a,b) ((a)>(b) ? (a):(b))




struct Subgrid{ 
	double *ex, *ey, *ez;              // electric-field unknowns
	double *hx, *hy, *hz;              // magnetic-field unknowns 
	double *cex, *cey, *cez;           // electric-field capacity update coefficient
	double *chx, *chy, *chz;           // magnetic-field capacity update coefficient
	double *sex, *sey, *sez;           // electric-field conductivity update coefficient
	double *aex, *aey, *aez;           // auxiliary time-average electric-field unknowns 
	double *qex, *qey, *qez;           // auxiliary electric-field coefficients
	double *dx, *dy, *dz;              // primary-grid steps
	int nx, ny, nz;                    // number of cells          	
	int bi[6];                         // bound indices in the main grid   
	int impl;                          // specifies which dimensions are solved implicitly   	
};   
typedef struct Subgrid Subgrid;  
                                  


// macros for accessing arrays and other attributes of the Subgrid struct s                       
#define Dx1_(I)        (*g).sg[s].dx[I]
#define Dy1_(J)        (*g).sg[s].dy[J]
#define Dz1_(K)        (*g).sg[s].dz[K]
#define Dx2_(I)        (0.5*(Dx1_(I)+Dx1_(I+1)))
#define Dy2_(J)        (0.5*(Dy1_(I)+Dy1_(J+1)))
#define Dz2_(K)        (0.5*(Dz1_(I)+Dz1_(K+1)))
#define Nx_            (*g).sg[s].nx
#define Ny_            (*g).sg[s].ny
#define Nz_            (*g).sg[s].nz
#define Bi_(I)         (*g).sg[s].bi[I]
#define Impl_          (*g).sg[s].impl
#define Ex_(I,J,K)     (*g).sg[s].ex[  ((I)*(Ny+1)+J)*(Nz+1)+K]  
#define Aex_(I,J,K)    (*g).sg[s].aex[ ((I)*(Ny+1)+J)*(Nz+1)+K]                            
#define Cex_(I,J,K)    (*g).sg[s].cex[ ((I)*(Ny-1)+J)*(Nz-1)+K]  
#define Sex_(I,J,K)    (*g).sg[s].sex[ ((I)*(Ny-1)+J)*(Nz-1)+K] 
#define Qex_(I,J,K,L)  (*g).sg[s].qex[(((I)*(Ny-1)+J)*(Nz-1)+K)*2+L]          
#define Ey_(I,J,K)     (*g).sg[s].ey[  ((I)*Ny+J)*(Nz+1)+K]
#define Aey_(I,J,K)    (*g).sg[s].aey[ ((I)*Ny+J)*(Nz+1)+K]
#define Cey_(I,J,K)    (*g).sg[s].cey[ ((I)*Ny+J)*(Nz-1)+K]
#define Sey_(I,J,K)    (*g).sg[s].sey[ ((I)*Ny+J)*(Nz-1)+K]
#define Qey_(I,J,K,L)  (*g).sg[s].qey[(((I)*Ny+J)*(Nz-1)+K)*2+L]
#define Ez_(I,J,K)     (*g).sg[s].ez[  ((I)*(Ny+1)+J)*Nz+K]
#define Aez_(I,J,K)    (*g).sg[s].aez[ ((I)*(Ny+1)+J)*Nz+K]
#define Cez_(I,J,K)    (*g).sg[s].cez[ ((I)*(Ny-1)+J)*Nz+K]
#define Sez_(I,J,K)    (*g).sg[s].sez[ ((I)*(Ny-1)+J)*Nz+K]
#define Qez_(I,J,K,L)  (*g).sg[s].qez[(((I)*(Ny-1)+J)*Nz+K)*2+L]
#define Hx_(I,J,K)     (*g).sg[s].hx[  ((I)*Ny+J)*Nz+K]
#define Chx_(I,J,K)    (*g).sg[s].chx[ ((I)*Ny+J)*Nz+K]
#define Hy_(I,J,K)     (*g).sg[s].hy[  ((I)*(Ny-1)+J)*Nz+K]
#define Chy_(I,J,K)    (*g).sg[s].chy[ ((I)*(Ny-1)+J)*Nz+K]
#define Hz_(I,J,K)     (*g).sg[s].hz[  ((I)*Ny+J)*(Nz-1)+K]
#define Chz_(I,J,K)    (*g).sg[s].chz[ ((I)*Ny+J)*(Nz-1)+K]


struct Schrodgrid {
	
	double *pr, *prv, *prvv;	   // real part of electron wavefunction at previous steps
	double *pi, *piv, *pivv;	   // imag part of electron wavefunction at previous steps
	double *lpx, *lpy, *lpz;	   // laplacian term update coefficient (positive)
	double *ppx, *ppy, *ppz;	   // potential term update coefficient (positive) -- for potential Vs 
	double *distance, *vs;		   // distance to center and static potential for the electron
	double *dxs, *dys, *dzs;	   // schrodinger grid steps -- may be the same as for Maxwell
	double *jx, *jy, *jz;		   // Quantum current unknowns
	int nx, ny, nz, nt, it; 	   //
	double dts;			   // timestep for Schrodinger, may be the same as for Maxwell

};
typedef struct Schrodgrid Schrodgrid;

// macros for accessing properties of the Schrodgrid
#define Dxs(I)	      (*g).srg[s].dxs[I]
#define Dys(J)	      (*g).srg[s].dys[J]
#define Dzs(K)	      (*g).srg[s].dzs[K]
#define Nxs			  (*g).srg[s].nx
#define Nys			  (*g).srg[s].ny
#define Nzs			  (*g).srg[s].nz
#define Nts			  (*g).srg[s].nt
#define Its			  (*g).srg[s].it
#define Dts			  (*g).srg[s].dts
#define Pr(I,J,K)     (*g).srg[s].pr[  ((I)*Nys+J)*(Nzs)+K]		
#define Prv(I,J,K)    (*g).srg[s].prv[ ((I)*Nys+J)*(Nzs)+K]		
#define Prvv(I,J,K)   (*g).srg[s].prvv[((I)*Nys+J)*(Nzs)+K]
#define Pi(I,J,K)     (*g).srg[s].pi[  ((I)*Nys+J)*(Nzs)+K]
#define Piv(I,J,K)    (*g).srg[s].piv[ ((I)*Nys+J)*(Nzs)+K]
#define Pivv(I,J,K)   (*g).srg[s].pivv[((I)*Nys+J)*(Nzs)+K]
#define Lpx(I)    	  (*g).srg[s].lpx[I]
#define Lpy(J)    	  (*g).srg[s].lpy[J]
#define Lpz(K)    	  (*g).srg[s].lpz[K]
#define Ppx(I)	      (*g).srg[s].ppx[I]
#define Ppy(J)  	  (*g).srg[s].ppy[J]
#define Ppz(K)  	  (*g).srg[s].ppz[K]
#define Vs(I,J,K)     (*g).srg[s].vs[((I)*Nys+J)*(Nzs)+K]
#define Dist(I,J,K)	  (*g).srg[s].distance[((I)*Nys+J)*(Nzs)+K]



struct Grid {
	double *ex, *ey, *ez;              // electric-field unknowns
	double *hx, *hy, *hz;              // magnetic-field unknowns 
	double *diffex, *diffey, *diffez;  // auxiliary time-difference electric-field unknowns 
	double *diffhx, *diffhy, *diffhz;  // auxiliary time-difference magnetic-field unknowns 
	double *cex, *cey, *cez;           // electric-field capacity update coefficient
	double *chx, *chy, *chz;           // magnetic-field capacity update coefficient
	double *sex, *sey, *sez;           // electric-field conductivity update coefficient
	double *dx, *dy, *dz;              // primary-grid steps
	double dtau;                       // Minkowski time step (c*dt)
	int nx, ny, nz, nt, it;            // number of cells, number of iterations, and current iteration  
	int bc[6];                         // boundary conditions          	
	int npml[6];                       // number of PML layers  
	int nsub, nsr;                     // number of subgrids and number of schrodinger grids  
	Subgrid *sg;                       // subgrids
	Schrodgrid *srg;		   // Schrodinger grids		   	    
}; 
typedef struct Grid Grid;           



          
// macros for accessing arrays and other attributes of the Grid struct g  
#define Dx1(I)        (*g).dx[I]
#define Dy1(J)        (*g).dy[J]
#define Dz1(K)        (*g).dz[K]
#define Dx2(I)        (0.5*(Dx1(I)+Dx1(I+1)))  
#define Dy2(J)        (0.5*(Dy1(J)+Dy1(J+1))) 
#define Dz2(K)        (0.5*(Dz1(K)+Dz1(K+1))) 
#define Dtau          (*g).dtau
#define Nx            (*g).nx
#define Ny            (*g).ny
#define Nz            (*g).nz
#define It            (*g).it
#define Nt            (*g).nt
#define BC(I)         (*g).bc[I] 
#define Npml(I)       (*g).npml[I] 
#define Nsub          (*g).nsub
#define Ex(I,J,K)     (*g).ex[    ((I)*(Ny+1)+J)*(Nz+1)+K]
#define DiffEx(I,J,K) (*g).diffex[((I)*(Ny-1)+J)*(Nz-1)+K]                              
#define Cex(I,J,K)    (*g).cex[   ((I)*(Ny-1)+J)*(Nz-1)+K]          
#define Sex(I,J,K)    (*g).sex[   ((I)*(Ny-1)+J)*(Nz-1)+K]            
#define Ey(I,J,K)     (*g).ey[    ((I)*Ny+J)*(Nz+1)+K]
#define DiffEy(I,J,K) (*g).diffey[((I)*Ny+J)*(Nz-1)+K]
#define Cey(I,J,K)    (*g).cey[   ((I)*Ny+J)*(Nz-1)+K]
#define Sey(I,J,K)    (*g).sey[   ((I)*Ny+J)*(Nz-1)+K]
#define Ez(I,J,K)     (*g).ez[    ((I)*(Ny+1)+J)*Nz+K]
#define DiffEz(I,J,K) (*g).diffez[((I)*(Ny-1)+J)*Nz+K]
#define Cez(I,J,K)    (*g).cez[   ((I)*(Ny-1)+J)*Nz+K]
#define Sez(I,J,K)    (*g).sez[   ((I)*(Ny-1)+J)*Nz+K]
#define Hx(I,J,K)     (*g).hx[    ((I)*Ny+J)*Nz+K]
#define DiffHx(I,J,K) (*g).diffhx[((I)*Ny+J)*Nz+K]
#define Chx(I,J,K)    (*g).chx[   ((I)*Ny+J)*Nz+K]
#define Hy(I,J,K)     (*g).hy[    ((I)*(Ny-1)+J)*Nz+K]
#define DiffHy(I,J,K) (*g).diffhy[((I)*(Ny-1)+J)*Nz+K]
#define Chy(I,J,K)    (*g).chy[   ((I)*(Ny-1)+J)*Nz+K]
#define Hz(I,J,K)     (*g).hz[    ((I)*Ny+J)*(Nz-1)+K]
#define DiffHz(I,J,K) (*g).diffhz[((I)*Ny+J)*(Nz-1)+K]
#define Chz(I,J,K)    (*g).chz[   ((I)*Ny+J)*(Nz-1)+K]

#define Nsr	      (*g).nsr 		// macro for accessing the number of schrodinger grids






// material mapping macros
// 'PyList_GET_ITEM' returns a borrowed reference. Hence, the reference count is okay as long as 'Py_DECREF'
// is applied to 'py_map' and 'py_lut', i.e. the objects that own the reference.
#define map(I,J,K)     PyInt_AS_LONG( PyList_GET_ITEM( py_map, (Py_ssize_t) ((I)*Ny+J)*Nz+K ) )
#define mur(I)         PyFloat_AS_DOUBLE( PyTuple_GET_ITEM( PyList_GET_ITEM( py_lut, (Py_ssize_t) I), 0) )
#define epsr(I)        PyFloat_AS_DOUBLE( PyTuple_GET_ITEM( PyList_GET_ITEM( py_lut, (Py_ssize_t) I), 1) )
#define sigma(I)       PyFloat_AS_DOUBLE( PyTuple_GET_ITEM( PyList_GET_ITEM( py_lut, (Py_ssize_t) I), 2) )
#define murx(I,J,K)    ((mur(map(I,J,K))+mur(map(I+1,J,K)))/2.0)
#define mury(I,J,K)    ((mur(map(I,J,K))+mur(map(I,J+1,K)))/2.0)
#define murz(I,J,K)    ((mur(map(I,J,K))+mur(map(I,J,K+1)))/2.0)
#define epsrx(I,J,K)   ((epsr(map(I,J,K))+epsr(map(I,J+1,K))+epsr(map(I,J,K+1))+epsr(map(I,J+1,K+1)))/4.0)
#define epsry(I,J,K)   ((epsr(map(I,J,K))+epsr(map(I+1,J,K))+epsr(map(I,J,K+1))+epsr(map(I+1,J,K+1)))/4.0)
#define epsrz(I,J,K)   ((epsr(map(I,J,K))+epsr(map(I+1,J,K))+epsr(map(I,J+1,K))+epsr(map(I+1,J+1,K)))/4.0)
#define sigmax(I,J,K)  (Dtau*Z0*(sigma(map(I,J,K))+sigma(map(I,J+1,K))+sigma(map(I,J,K+1))+sigma(map(I,J+1,K+1)))/4.0)
#define sigmay(I,J,K)  (Dtau*Z0*(sigma(map(I,J,K))+sigma(map(I+1,J,K))+sigma(map(I,J,K+1))+sigma(map(I+1,J,K+1)))/4.0)
#define sigmaz(I,J,K)  (Dtau*Z0*(sigma(map(I,J,K))+sigma(map(I+1,J,K))+sigma(map(I,J+1,K))+sigma(map(I+1,J+1,K)))/4.0)



// function prototypes
PyObject* run_fdtd(PyObject *ipt); 

void initGrid(Grid *g, PyObject *ipt); 
void updateE(Grid *g);
void updateH(Grid *g); 
void updateDiffE(Grid *g);
void updateDiffH(Grid *g);   
void freeMemoryGrid(Grid *g); 

void initSubgrids(Grid *g, PyObject *ipt);
void updateSubgrids(Grid *g);
void fine2coarse(Grid *g);
void coarse2fine(Grid *g); 
void freeMemorySubgrids(Grid *g);

void initSource(Grid *g, PyObject *ipt); 
void updateDiffEcurrent(Grid *g);
void updateEhard(Grid *g);
void updateHhard(Grid *g);
void freeMemorySource(void); 

void initSensor(Grid *g, PyObject *ipt);
void updateSensor(Grid *g);
PyObject* sensorC2Py(Grid *g);

void initPml(Grid *g, PyObject *ipt);	
void updateDiffEpml(Grid *g);
void updateDiffHpml(Grid *g);
void freeMemoryPml(void);

void initAdhie1(Grid *g, PyObject *ipt);
void initAdhie2(Grid *g);
void updateDiffEadhie(Grid *g);
void updateDiffHadhie(Grid *g);
void freeMemoryAdhie(void); 

void initVisualization(Grid *g, PyObject *ipt); 
void visualize(Grid *g); 
void freeMemoryVisualization(void);

void initSchrodgrid(Grid *g, PyObject *ipt);	// new function prototypes start here
void freeMemorySchrodgrid(Grid *g);
void electronSnapshot(Grid *g);
void updateElectron(Grid *g);
void electronTimeevolution(Grid *g);





#endif   // matches #ifndef _FDTD3D_H in the beginning of the header file



/* COMMENTS: 
(1) The requirement of memory compactness, i.e. all array elements need to 
      be stored as close as possible to each other in memory, in combination 
      with dynamic memory allocation (since Nx, Ny and Nz are not known at 
      compile time) makes this the most efficient way to implement the 3d
      arrays. As explained in Alain's "Jumping into C++", pointers to 
      pointers can be arbitrarily far in memory when allocated dynamically.  
(2) Nodes that are adjacent to each other in the z-direction are also 
      adjacent to each other in memory (row-major vectorization).
      Consequently, the inner loops should iterate over the z-index.
      Row-major ordering is the default in C/C++ and numpy.  
(3) The order of freeing memory must be the reverse order of allocating
      memory to avoid memory leaks. First traverse to the child memory 
      location and start freeing from there, traversing back to the parent 
      node. Also, note that 'free(malloc(null))' is valid.
(4) Arrays are automatically dereferenced (the value is accessed instead of
      the address) by the "[ ]" operator. No '*' is needed. 
(5) Cython does not recognize "->".
(6) Convention used here: Macros that access Grid and Subgrid elements start 
      with a capital letter. Subgrid macros have a trailing underscore. 
      Mind that the leading underscore is reserved by the compiler. 
(7) Be very careful when using g and s to denote other things than the 
      main grid object or the subgrid index.
(8) Perhaps, single instead of double precision could be used since FDTD
      is quite robust against round-off errors. However, this limits the 
      the multiscale capabilities of the solver. 
*/

 

