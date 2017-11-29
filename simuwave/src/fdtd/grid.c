
#include "fdtd3d.h"



void initGrid(Grid *g, PyObject *ipt){	

	int i, j, k; 

	// number of iterations
	PyObject *py_Nt = PyObject_GetAttrString( ipt, "Nt" );
	Nt = PyInt_AS_LONG(py_Nt);
	Py_DECREF(py_Nt);
	
	// time step
	PyObject *py_mg = PyObject_GetAttrString( ipt, "maingrid" ); 	
	PyObject *py_dtau = PyObject_GetAttrString( py_mg, "dtau" ); // dtau = c*dt
	Dtau = 1e-3 * PyFloat_AS_DOUBLE(py_dtau);  // (m) -> miliseconds				
	Py_DECREF(py_dtau);

	// number of cells	
	PyObject *py_steps = PyObject_GetAttrString( py_mg, "steps" );   // new reference => Py_DECREF
	PyObject *py_dx = PyList_GET_ITEM(py_steps, 0);                 // borrowed reference => no Py_DECREF
	PyObject *py_dy = PyList_GET_ITEM(py_steps, 1); 
	PyObject *py_dz = PyList_GET_ITEM(py_steps, 2);
	Nx = PyList_Size(py_dx);	
	Ny = PyList_Size(py_dy); 
	Nz = PyList_Size(py_dz); 	

	// cell dimensions
	my_alloc( (*g).dx, Nx , double );  
	my_alloc( (*g).dy, Ny , double ); 	
	my_alloc( (*g).dz, Nz , double ); 
	for(i=0; i<Nx; i++){ Dx1(i) = 1e-3 * PyFloat_AS_DOUBLE( PyList_GET_ITEM( py_dx, (Py_ssize_t) i ) ); }  // (m)
	for(j=0; j<Ny; j++){ Dy1(j) = 1e-3 * PyFloat_AS_DOUBLE( PyList_GET_ITEM( py_dy, (Py_ssize_t) j ) ); }	
	for(k=0; k<Nz; k++){ Dz1(k) = 1e-3 * PyFloat_AS_DOUBLE( PyList_GET_ITEM( py_dz, (Py_ssize_t) k ) ); }		
	Py_DECREF(py_steps);

	// medium (needed for the macros defined in the header file)
	PyObject *py_map = PyObject_GetAttrString(py_mg, "matmap"); 	
	PyObject *py_lut = PyObject_GetAttrString(py_mg, "matlut"); 

	// dynamically allocate memory for the field variables and update coefficients (intialized to zero)
	my_alloc( (*g).ex     , Nx*(Ny+1)*(Nz+1) , double );  
	my_alloc( (*g).diffex , Nx*(Ny-1)*(Nz-1) , double ); 	
	my_alloc( (*g).cex    , Nx*(Ny-1)*(Nz-1) , double ); 
	my_alloc( (*g).sex    , Nx*(Ny-1)*(Nz-1) , double ); 
        my_alloc( (*g).ey     , (Nx+1)*Ny*(Nz+1) , double ); 
	my_alloc( (*g).diffey , (Nx-1)*Ny*(Nz-1) , double ); 	
	my_alloc( (*g).cey    , (Nx-1)*Ny*(Nz-1) , double );
	my_alloc( (*g).sey    , (Nx-1)*Ny*(Nz-1) , double );
	my_alloc( (*g).ez     , (Nx+1)*(Ny+1)*Nz , double ); 
	my_alloc( (*g).diffez , (Nx-1)*(Ny-1)*Nz , double ); 
	my_alloc( (*g).cez    , (Nx-1)*(Ny-1)*Nz , double ); 
	my_alloc( (*g).sez    , (Nx-1)*(Ny-1)*Nz , double ); 
 	my_alloc( (*g).hx     , (Nx-1)*Ny*Nz     , double ); 
	my_alloc( (*g).diffhx , (Nx-1)*Ny*Nz     , double ); 
	my_alloc( (*g).chx    , (Nx-1)*Ny*Nz     , double ); 
	my_alloc( (*g).hy     , Nx*(Ny-1)*Nz     , double );
	my_alloc( (*g).diffhy , Nx*(Ny-1)*Nz     , double );
	my_alloc( (*g).chy    , Nx*(Ny-1)*Nz     , double );
	my_alloc( (*g).hz     , Nx*Ny*(Nz-1)     , double );
	my_alloc( (*g).diffhz , Nx*Ny*(Nz-1)     , double );
	my_alloc( (*g).chz    , Nx*Ny*(Nz-1)     , double );

	
	


	// set electric-field update coefficients 	
	for(i=0; i<Nx; i++){  
	  for(j=0; j<Ny-1; j++){
	    for(k=0; k<Nz-1; k++){  	
	      Cex(i,j,k) =  Dtau * Dx1(i) / ( Dy2(j) * Dz2(k) * ( epsrx(i,j,k) + sigmax(i,j,k) ) );	
	      Sex(i,j,k) = -sigmax(i,j,k) / ( epsrx(i,j,k) + sigmax(i,j,k) );	 		      
	}}}	
 
	for(i=0; i<Nx-1; i++){  
	  for(j=0; j<Ny; j++){
	    for(k=0; k<Nz-1; k++){           	
	      Cey(i,j,k) =  Dtau * Dy1(j) / ( Dx2(i) * Dz2(k) * (epsry(i,j,k)+sigmay(i,j,k)) );		
	      Sey(i,j,k) = -sigmay(i,j,k) / ( epsry(i,j,k) + sigmay(i,j,k) );         
	}}}	

	for(i=0; i<Nx-1; i++){  
	  for(j=0; j<Ny-1; j++){
	    for(k=0; k<Nz; k++){ 
	      Cez(i,j,k) =  Dtau * Dz1(k) / ( Dx2(i) * Dy2(j) * (epsrz(i,j,k)+sigmaz(i,j,k)) ); 
	      Sez(i,j,k) = -sigmaz(i,j,k) / ( epsrz(i,j,k) + sigmaz(i,j,k) );         			      
	}}}	

	// set magnetic-field update coefficients
	for(i=0; i<Nx-1; i++){  
	  for(j=0; j<Ny; j++){
	    for(k=0; k<Nz; k++){
	      Chx(i,j,k) = -Dtau * Dx2(i) / ( Dy1(j) * Dz1(k) * murx(i,j,k) );   
        }}}		

	for(i=0; i<Nx; i++){  
	  for(j=0; j<Ny-1; j++){
	    for(k=0; k<Nz; k++){	         
	      Chy(i,j,k) = -Dtau * Dy2(j) / ( Dx1(i) * Dz1(k) * mury(i,j,k) );   
	}}}	

	for(i=0; i<Nx; i++){  
	  for(j=0; j<Ny; j++){
	    for(k=0; k<Nz-1; k++){	         
              Chz(i,j,k) = -Dtau * Dz2(k) / ( Dx1(i) * Dy1(j) * murz(i,j,k) );    
	}}}	

	// decrease reference counts	
	Py_DECREF(py_map); 
        Py_DECREF(py_lut); 
	Py_DECREF(py_mg);

}



void updateDiffE(Grid *g){
	
	int i, j, k; 

	for(i=0; i<Nx; i++){  
	  for(j=0; j<Ny-1; j++){
	    for(k=0; k<Nz-1; k++){              		
	      DiffEx(i,j,k) = Cex(i,j,k) * ((Hz(i,j+1,k)-Hz(i,j,k)) - (Hy(i,j,k+1)-Hy(i,j,k))) + Sex(i,j,k) * Ex(i,j+1,k+1); 
	}}}

	for(i=0; i<Nx-1; i++){  
	  for(j=0; j<Ny; j++){
	    for(k=0; k<Nz-1; k++){ 		
	      DiffEy(i,j,k) = Cey(i,j,k) * ((Hx(i,j,k+1)-Hx(i,j,k)) - (Hz(i+1,j,k)-Hz(i,j,k))) + Sey(i,j,k) * Ey(i+1,j,k+1);
	}}}
	
	for(i=0; i<Nx-1; i++){  
	  for(j=0; j<Ny-1; j++){
	    for(k=0; k<Nz; k++){ 		
	      DiffEz(i,j,k) = Cez(i,j,k) * ((Hy(i+1,j,k)-Hy(i,j,k)) - (Hx(i,j+1,k)-Hx(i,j,k))) + Sez(i,j,k) * Ez(i+1,j+1,k); 
	}}}

}



void updateDiffH(Grid *g){
	
	int i, j, k; 

	for(i=0; i<Nx-1; i++){  
	  for(j=0; j<Ny; j++){
	    for(k=0; k<Nz; k++){		
	      DiffHx(i,j,k) = Chx(i,j,k) * ((Ez(i+1,j+1,k)-Ez(i+1,j,k)) - (Ey(i+1,j,k+1)-Ey(i+1,j,k)));	
	}}}

	for(i=0; i<Nx; i++){  
	  for(j=0; j<Ny-1; j++){
	    for(k=0; k<Nz; k++){				
	      DiffHy(i,j,k) = Chy(i,j,k) * ((Ex(i,j+1,k+1)-Ex(i,j+1,k)) - (Ez(i+1,j+1,k)-Ez(i,j+1,k))); 
	}}}

	for(i=0; i<Nx; i++){  
	  for(j=0; j<Ny; j++){
	    for(k=0; k<Nz-1; k++){ 				
	      DiffHz(i,j,k) = Chz(i,j,k) * ((Ey(i+1,j,k+1)-Ey(i,j,k+1)) - (Ex(i,j+1,k+1)-Ex(i,j,k+1)));
	}}}

}



void updateE(Grid *g){
	
	int i, j, k; 

	for(i=0; i<Nx; i++){  
	  for(j=0; j<Ny-1; j++){
	    for(k=0; k<Nz-1; k++){              		
	      Ex(i,j+1,k+1) += DiffEx(i,j,k); 
	}}}

	for(i=0; i<Nx-1; i++){  
	  for(j=0; j<Ny; j++){
	    for(k=0; k<Nz-1; k++){ 		
	      Ey(i+1,j,k+1) += DiffEy(i,j,k);
	}}}
	
	for(i=0; i<Nx-1; i++){  
	  for(j=0; j<Ny-1; j++){
	    for(k=0; k<Nz; k++){ 		
	      Ez(i+1,j+1,k) += DiffEz(i,j,k); 
	}}}

}



void updateH(Grid *g){
	
	int i, j, k; 

	for(i=0; i<Nx-1; i++){  
	  for(j=0; j<Ny; j++){
	    for(k=0; k<Nz; k++){		
	      Hx(i,j,k) += DiffHx(i,j,k);	
	}}}

	for(i=0; i<Nx; i++){  
	  for(j=0; j<Ny-1; j++){
	    for(k=0; k<Nz; k++){				
	      Hy(i,j,k) += DiffHy(i,j,k); 
	}}}

	for(i=0; i<Nx; i++){  
	  for(j=0; j<Ny; j++){
	    for(k=0; k<Nz-1; k++){ 				
	      Hz(i,j,k) += DiffHz(i,j,k);
	}}}

}



void freeMemoryGrid(Grid *g){     
                            
	                free((*g).chz); free((*g).diffhz); free((*g).hz);	
	                free((*g).chy); free((*g).diffhy); free((*g).hy); 		
	                free((*g).chx); free((*g).diffhx); free((*g).hx);     	
	free((*g).sez); free((*g).cez); free((*g).diffez); free((*g).ez);
	free((*g).sey); free((*g).cey); free((*g).diffey); free((*g).ey);	
	free((*g).sex); free((*g).cex); free((*g).diffex); free((*g).ex); 	
	free((*g).dz); free((*g).dy); free((*g).dx);    
	free(g); 

}



/* COMMENTS:
(1) Use normalized fields (E_num=E_phys/Z0) and time step (Dtau=c0*dt).
(2) Rescale the normalized electric and magnetic fields by the local step size as to actually yield 
      voltages and currents. This reduces the flops, the memory queries and the memory usage. 
(3) Use (zero) dummy variables for the tangential electric fields at the exterior boundaries. 
(4) The auxiliary difference unknowns are defined as DiffEx(n) = Ex((n+1)*dt)-Ex(n*dt)
(5) We use the time-forward discretization of the conduction current because this reduces the ADHIE error
     at negligible accuracy cost. 
(6) Multiplications are faster than divisions, so avoid divisions in the update equations.
*/




