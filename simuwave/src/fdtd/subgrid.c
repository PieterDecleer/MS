
#include "subgrid.h" 


void initSubgrids(Grid *g, PyObject *ipt){		

	// three spatial indices and subgrid index
	int i, j, k, s;  

	// number of subgrids
	PyObject *py_sg = PyObject_GetAttrString( ipt, "subgrid" ); 	
	Nsub = PyList_Size(py_sg);
	

	for(s=0; s<Nsub; s++){	
	 
	  PyObject *py_sg_s = PyList_GET_ITEM(py_sg,s);	

	  // number of cells
	  PyObject *py_steps = PyObject_GetAttrString(py_sg, "steps");
	  PyObject *py_Dx_ = PyList_GET_ITEM(py_steps, 0); 
	  PyObject *py_Dy_ = PyList_GET_ITEM(py_steps, 1); 
	  PyObject *py_Dz_ = PyList_GET_ITEM(py_steps, 2); 
	  Nx_ = PyList_Size(py_Dx_);	
	  Ny_ = PyList_Size(py_Dy_); 
	  Nz_ = PyList_Size(py_Dz_);	

	  // cell dimensions
	  my_alloc( (*g).sg[s].dx, Nx_ , double );  
	  my_alloc( (*g).sg[s].dy, Ny_ , double ); 
	  my_alloc( (*g).sg[s].dz, Nz_ , double );  
	  for(i=0; i<Nx_; i++){ Dx1_(i) = 1e-3 * PyFloat_AS_DOUBLE( PyList_GET_ITEM( py_Dx_, (Py_ssize_t) i ) ); }   // (m)
	  for(j=0; j<Ny_; j++){ Dy1_(j) = 1e-3 * PyFloat_AS_DOUBLE( PyList_GET_ITEM( py_Dy_, (Py_ssize_t) j ) ); } 
      	  for(k=0; k<Nz_; k++){ Dz1_(k) = 1e-3 * PyFloat_AS_DOUBLE( PyList_GET_ITEM( py_Dz_, (Py_ssize_t) k ) ); }   
	  Py_DECREF(py_steps);   

	  // implicitization direction
	  PyObject *py_impl = PyObject_GetAttrString( py_sg, "impl" );		 
	  Impl_=0; for (i=0; i<3; i++){ Impl_ += pow(2,i) * PyInt_AS_LONG( PyList_GET_ITEM( py_impl , (Py_ssize_t) i ) ); }
          Py_DECREF(py_impl);

	  // subgrid-bounds indices
	  PyObject *py_bi = PyObject_GetAttrString( py_sg, "bi" );	
	  for (i=0; i<6; i++){ Bi_(i) = PyInt_AS_LONG(PyList_GET_ITEM( py_bi , (Py_ssize_t) i ) ); }
	  Py_DECREF(py_bi);

	  // medium (needed for the macros defined in the header file)
	  PyObject *py_map = PyObject_GetAttrString( py_sg, "matmap" ); 	
	  PyObject *py_lut = PyObject_GetAttrString( py_sg, "matlut" ); 

	  // dynamically allocate memory for the field variables and update coefficients (intialized to zero)
	  my_alloc( (*g).sg[s].ex  , Nx_*(Ny_+1)*(Nz_+1) , double );  
	  my_alloc( (*g).sg[s].cex , Nx_*(Ny_-1)*(Nz_-1) , double );
	  my_alloc( (*g).sg[s].sex , Nx_*(Ny_-1)*(Nz_-1) , double );
	  my_alloc( (*g).sg[s].ey  , (Nx_+1)*Ny_*(Nz_+1) , double ); 
	  my_alloc( (*g).sg[s].cey , (Nx_-1)*Ny_*(Nz_-1) , double );
	  my_alloc( (*g).sg[s].sey , (Nx_-1)*Ny_*(Nz_-1) , double );
	  my_alloc( (*g).sg[s].ez  , (Nx_+1)*(Ny_+1)*Nz_ , double ); 
	  my_alloc( (*g).sg[s].cez , (Nx_-1)*(Ny_-1)*Nz_ , double ); 
	  my_alloc( (*g).sg[s].sez , (Nx_-1)*(Ny_-1)*Nz_ , double ); 
	  my_alloc( (*g).sg[s].hx  , (Nx_-1)*Ny_*Nz_     , double );  
	  my_alloc( (*g).sg[s].chx , (Nx_-1)*Ny_*Nz_     , double ); 
	  my_alloc( (*g).sg[s].hy  , Nx_*(Ny_-1)*Nz_     , double );
	  my_alloc( (*g).sg[s].chy , Nx_*(Ny_-1)*Nz_     , double );
	  my_alloc( (*g).sg[s].hz  , Nx_*Ny_*(Nz_-1)     , double );
	  my_alloc( (*g).sg[s].chz , Nx_*Ny_*(Nz_-1)     , double );


	  // set standard electric-field update coefficients 	
	  for(i=0; i<Nx_; i++){  
	   for(j=0; j<Ny_-1; j++){
	    for(k=0; k<Nz_-1; k++){  
             temp = ;        	
	     Cex_(i,j,k) = Dtau * Dx1_(i) / ( Dy2_(j) * Dz2_(k) * ( epsrx(i,j,k) + sigmax(i,j,k)/2.0 ) ); 
	     Sex_(i,j,k) = ( epsrx(i,j,k) - sigmax(i,j,k)/2.0 ) / ( epsrx(i,j,k) + sigmax(i,j,k)/2.0 );  		      
	   }}}	

	  for(i=0; i<Nx_-1; i++){  
	   for(j=0; j<Ny_; j++){
	    for(k=0; k<Nz_-1; k++){       	
             Cey_(i,j,k) = Dtau * Dy1_(j) / ( Dx2_(i) * Dz2_(k) * ( epsry(i,j,k) + sigmay(i,j,k)/2.0 ) );		
	     Sey_(i,j,k) = ( epsry(i,j,k) - sigmay(i,j,k)/2.0 ) / ( epsry(i,j,k) + sigmay(i,j,k)/2.0 );         
	  }}}	

	  for(i=0; i<Nx_-1; i++){  
	   for(j=0; j<Ny_-1; j++){
	    for(k=0; k<Nz_; k++){ 
	     Cez_(i,j,k) = Dtau * Dz1_(k) / ( Dx2_(i) * Dy2_(j) * ( epsrz(i,j,k) + sigmaz(i,j,k)/2.0 ) ); 
	     Sez_(i,j,k) = ( epsrz(i,j,k) - sigmaz(i,j,k)/2.0 ) / ( epsrz(i,j,k) + sigmaz(i,j,k)/2.0 );
	  }}}	

	  // set standard magnetic-field update coefficients
	  for(i=0; i<Nx_-1; i++){  
	   for(j=0; j<Ny_; j++){
	    for(k=0; k<Nz_; k++){             
	     Chx_(i,j,k) = -Dtau * Dx2_(i) / ( Dy1_(j) * Dz1_(k) * murx(i,j,k) );     
	  }}}		

	  for(i=0; i<Nx_; i++){  
	   for(j=0; j<Ny_-1; j++){
	    for(k=0; k<Nz_; k++){
	     Chy_(i,j,k) = -Dtau * Dy2_(j) / ( Dx1_(i) * Dz1_(k) * mury(i,j,k) );  
	  }}}	

	  for(i=0; i<Nx_; i++){  
	   for(j=0; j<Ny_; j++){
	    for(k=0; k<Nz_-1; k++){	         
	     Chz_(i,j,k) = -Dtau * Dz2_(k) / ( Dx1_(i) * Dy1_(j) * murz(i,j,k) );    
	  }}}	


	  // decrease reference counts	
	  Py_DECREF(py_map); 
	  Py_DECREF(py_lut); 
	  Py_DECREF(py_sg_s);
	  

          // preprocessing
	  switch (Impl_){

            case 1: // z
              break;

            case 2: // y
	
	      my_alloc( (*g).sg[s].aex , Nx_*(Ny_+1)*(Nz_+1)   , double );
	      my_alloc( (*g).sg[s].aez , (Nx_+1)*(Ny_+1)*Nz_   , double ); 		
	      my_alloc( (*g).sg[s].qex , Nx_*(Ny_-1)*(Nz_-1)*3 , double );  
	      my_alloc( (*g).sg[s].qez , (Nx_-1)*(Ny_-1)*Nz_*3 , double );	

	      for(i=0; i<Nx_; i++){  
	       for(j=0; j<Ny_-1; j++){
	        for(k=0; k<Nz_-1; k++){  		      
	         Qex_(i,j,k,0) = Cex_(i,j,k)*Chz_(i,j  ,k);                            // subdiagonal (negative)
	         Qex_(i,j,k,2) = Cex_(i,j,k)*Chz_(i,j+1,k);                            // superdiagonal (negative)   
	         Qex_(i,j,k,1) = 1.0-Qex_(i,j,k,0)-Qex_(i,j,k,2);                      // main diagonal (positive)
                 Qex_(i,j,k,1) = 1.0/(Qex_(i,j,k,1)-Qex_(i,j,k,0)*Qex_(i,j-1,k,2));    // decomposition (Thomas algorithm) 
	         Qex_(i,j,k,2) *= Qex_(i,j,k,1);
	      }}}
	      
	      for(i=0; i<Nx_-1; i++){  
	       for(j=0; j<Ny_-1; j++){
	        for(k=0; k<Nz_; k++){   		      
	         Qez_(i,j,k,0) = Cez_(i,j,k)*Chx_(i,j  ,k);                    
	         Qez_(i,j,k,2) = Cez_(i,j,k)*Chx_(i,j+1,k);                    
	         Qez_(i,j,k,1) = 1.0-Qez_(i,j,k,0)-Qez_(i,j,k,2);                       
                 Qez_(i,j,k,1) = 1.0/(Qez_(i,j,k,1)-Qez_(i,j,k,0)*Qez_(i,j-1,k,2));     
	         Qez_(i,j,k,2) *= Qez_(i,j,k,1);
	      }}}

	      break;

            case 3: // yz
	      break;

	    case 4: // x 
              break;
                                                         
	    case 5: // xz
	      break;

            case 6: // xy
	      break;

	    case 7: // xyz  
	      break;
                                          
	  }
	
	}   


	Py_DECREF(py_sg);                 
}



void updateSubgrids(Grid *g){       
	
	int i, j, k, s;             

	for(s=0; s<Nsub; s++){ 	
		  
  	  switch (Impl){
           
	    case 1: // z
              break;


            case 2: // y


	      // Ex 
	      for(i=0; i<Nx_; i++){  
	       for(j=0; j<Ny_-1; j++){
	        for(k=0; k<Nz_-1; k++){         	
		 _Aex(i,j+1,k+1) =  Cex_(i,j,k)*(Hz_(i,j+1,k)-Hz_(i,j,k)) - Cex_(i,j,k)*(Hy_(i,j,k+1)-Hy_(i,j,k)) + (1.0+Sex_(i,j,k))*Ex_(i,j+1,k+1)
				      + 0.5*Cex_(i,j,k)*( Chz_(i,j+1,k)*(Ey_(i+1,j+1,k+1)-Ey_(i,j+1,k+1)) - Chz_(i,j,k)*(Ey_(i+1,j,k+1)-Ey_(i,j,k+1)) );		      
	      }}}	
	      for (i=0; i<Nx; i++){  
	       for (k=0; k<Nz-1; k++){
	        _Aex(i,1,k+1) *= Qex_(i,0,k,1);
	        for (j=1    ; j<Ny_; j++){ _Aex(i,j+1,k+1) = (_Aex(i,j+1,k+1)-Qex_(i,j,k,0)*_Aex(i,j,k+1))*Qex_(i,j,k,1); }   // forward sweep                     
	        for (j=Ny_-2; j>=0 ; j--){ _Aex(i,j+1,k+1) -= Qex_(i,j,k,2)*_Aex(i,j+2,k+1); }                                // backward sweep 	
	      }}		     	  
	      for(i=0; i<Nx_; i++){  
	       for(j=0; j<Ny_-1; j++){
	        for(k=0; k<Nz_-1; k++){         	
		 Ex_(i,j+1,k+1) = _Aex(i,j+1,k+1)-Ex_(i,j+1,k+1);		      
	      }}}		


	      // Ez
	      for(i=0; i<Nx_-1; i++){  
	       for(j=0; j<Ny_-1; j++){
	        for(k=0; k<Nz_; k++){         	
		 _Aez(i+1,j+1,k) = Cez_(i,j,k)*(Hy_(i+1,j,k)-Hy_(i,j,k)) - Cez_(i,j,k)*(Hx_(i,j+1,k)-Hx_(i,j,k)) + (1.0+Sez_(i,j,k))*Ez_(i+1,j+1,k)      
		                     + 0.5*Cez_(i,j,k)*( Chx_(i,j+1,k)*(Ey_(i+1,j+1,k+1)-Ey_(i+1,j+1,k)) - Chx_(i,j,k)*(Ey_(i+1,j,k+1)-Ey_(i+1,j,k)) ); 			      
	      }}}		
	      for(i=0; i<Nx_-1; i++){  
	       for (k=0; k<Nz_; k++){
	        _Aez(i+1,1,k,3) *= Qez_(i,0,k,1);
	        for (j=1    ; j<Ny_; j++){ _Aez(i+1,j+1,k) = (_Aez(i+1,j+1,k)-Qez_(i,j,k,0)*_Aez(i+1,j,k))*Qez_(i,j,k,1); }   // forward sweep                    
	        for (j=Ny_-2; j>=0 ; j--){ _Aez(i+1,j+1,k) -= Qez_(i,j,k,2)*_Aez(i+1,j+2,k); }                                // backward sweep 	
	      }} 
	      for(i=0; i<Nx_-1; i++){  
	       for(j=0; j<Ny_-1; j++){
	        for(k=0; k<Nz_; k++){           	
		 Ez_(i+1,j+1,k) = _Aez(i+1,j+1,k)-Ez_(i+1,j+1,k);		      
	      }}}
	 	

	      // Hx
	      for(i=0; i<Nx_-1; i++){  
	       for(j=0; j<Ny_; j++){
	        for(k=0; k<Nz_; k++){		
	         Hx_(i,j,k) += Chx_(i,j,k) * (0.5*(_Aez(i+1,j+1,k)-_Aez(i+1,j,k)) - (Ey_(i+1,j,k+1)-Ey_(i+1,j,k)));	
 	      }}}


	      // Hy	
	      for(i=0; i<Nx_; i++){  
	       for(j=0; j<Ny_-1; j++){
	        for(k=0; k<Nz_; k++){				
	         Hy_(i,j,k) += Chy_(i,j,k) * ((Ex_(i,j+1,k+1)-Ex_(i,j+1,k)) - (Ez_(i+1,j+1,k)-Ez_(i,j+1,k))); 
	      }}}


              // Hz		
 	      for(i=0; i<Nx_; i++){  
	       for(j=0; j<Ny_; j++){
	        for(k=0; k<Nz_-1; k++){ 				
	         Hz_(i,j,k) += Chz_(i,j,k) * ((Ey_(i+1,j,k+1)-Ey_(i,j,k+1)) - 0.5*(_Aex(i,j+1,k+1)-_Aex(i,j,k+1)));
	      }}}

	      // Ey
	      for(i=0; i<Nx_-1; i++){  
	       for(j=0; j<Ny_; j++){
	        for(k=0; k<Nz_-1; k++){ 		
	         Ey_(i+1,j,k+1) = Cey_(i,j,k) * ((Hx_(i,j,k+1)-Hx_(i,j,k)) - (Hz_(i+1,j,k)-Hz_(i,j,k))) + Sey_(i,j,k) * Ey_(i+1,j,k+1);
	      }}}

	      break;



            case 3: // yz
   	      break;


	    case 4: // x 
              break;

                                                         
	    case 5: // xz
	      break;


            case 6: // xy
	      break;


	    case 7: // xyz  

                                          
	  }




void fine2coarse(Grid *g){         // complete the coarse updates of the auxiliary time-difference field unknowns

	int i, j, k, i_, j_, k_, s; 

	for(s=0; s<Nsub; s++){ 
	   
	  // x-faces  
	  for(j=0; j<Ny; j++){
           for(j_=...; j_<...; j_++){
	    for(k=0; k<Nz-1; k++){ 
	     for(k_=...; k_<...; k_++){				
	      DiffEy(Bi_(0)-1,j,k) += Cey(Bi_(0)-1,j,k) * ... * Hz_(0    ,j_,k_);
              DiffEy(Bi_(3)  ,j,k) -= Cey(Bi_(3)  ,j,k) * ... * Hz_(Nx_-1,j_,k_);
	  }}}}
	  for(j=0; j<Ny-1; j++){
           for(j_=...; j_<...; j_++){
	    for(k=0; k<Nz; k++){ 	
	     for(k_=...; k_<...; k_++){
	      DiffEz(Bi_(0)-1,j,k) -= Cez(Bi_(0)-1,j,k) * ... * Hy_(0    ,j_,k_); 
              DiffEz(Bi_(3)  ,j,k) += Cez(Bi_(3)  ,j,k) * ... * Hy_(Nx_-1,j_,k_); 
	  }}}}

          // y-faces	  
          for(i=0; i<Nx; i++){ 
           for(i_=...; i_<...; i_++){ 
	    for(k=0; k<Nz-1; k++){  
             for(k_=...; k_<...; k_++){             		
	      DiffEx(i,Bi_(1)-1,k) -= Cex(i,Bi_(1)-1,k) * ... * Hz_(i_,0    ,k_); 
              DiffEx(i,Bi_(4)  ,k) += Cex(i,Bi_(4)  ,k) * ... * Hz_(i_,Ny_-1,k_); 
	  }}}}
	  for(i=0; i<Nx-1; i++){ 
           for(i_=...; i_<...; i_++){ 
	    for(k=0; k<Nz; k++){ 
             for(k_=...; k_<...; k_++){		
	      DiffEz(i,Bi_(1)-1,k) += Cez(i,j,Bi_(1)-1) * ... * Hx_(i_,0    ,k_); 
              DiffEz(i,Bi_(4)  ,k) -= Cez(i,j,Bi_(4)  ) * ... * Hx_(i_,Ny_-1,k_); 
	  }}}}

	  // z-faces
          for(i=0; i<Nx; i++){  
           for(i_=...; i_<...; i_++){  
	    for(j=0; j<Ny-1; j++){ 
             for(j_=...; j_<...; j_++){             		
	      DiffEx(i,j,Bi_(2)-1) += Cex(i,j,Bi_(2)-1) * ... * Hy_(i_,j_,0    ); 
              DiffEx(i,j,Bi_(5)  ) -= Cex(i,j,Bi_(5)  ) * ... * Hy_(i_,j_,Nz_-1); 
	  }}}}
          for(i=0; i<Nx-1; i++){ 
           for(i_=...; i_<...; i_++){ 
            for(j=0; j<Ny; j++){
             for(j_=...; j_<...; j_++){ 		
	      DiffEy(i,j,Bi_(2)-1) -= Cey(i,j,Bi_(2)-1) * ... * Hx_(i_,j_,0    );
              DiffEy(i,j,Bi_(5)  ) += Cey(i,j,Bi_(5)  ) * ... * Hx_(i_,j_,Nz_-1);
	  }}}}	

	}
}



void coarse2fine(Grid *g){      // subgrid boundary condition

	int i, j, k, i_, j_, k_, s; 

	for(s=0; s<Nsub; s++){
 	  			
	  // x-faces  
	  for(j_=0; j_<Ny_; j_++){
	   for(k_=0; k_<Nz_-1; k_++){ 		
	    j = ... ;
	    k = ... ;
	    Ey_(0  ,j_,k_+1) = Dy1_(j_)/Dy1(j)*Ey(Bi_(0)  ,j,k+1);
            Ey_(Nx_,j_,k_+1) = Dy1_(j_)/Dy1(j)*Ey(Bi_(3)+1,j,k+1);
	  }}  
	  for(j_=0; j_<Ny-1; j_++){
	   for(k_=0; k_<Nz; k_++){ 		
	    j = ... ;
	    k = ... ;
	    Ez_(0  ,j_+1,k_) = Dz1_(k_)/Dz1(k)*Ez(Bi_(0)  ,j+1,k); 
	    Ez_(Nx_,j_+1,k_) = Dz1_(k_)/Dz1(k)*Ez(Bi_(3)+1,j+1,k);
	  }}

          // y-faces	  
          for(i_=0; i_<Nx; i_++){  
	   for(k_=0; k_<Nz-1; k_++){ 
	    i = ... ;
	    k = ... ;             		
	    Ex_(i_,0  ,k_+1) = Dx_(i_)/Dx1(i)*Ex(i,Bi_(1)  ,k+1); 
            Ex_(i_,Ny_,k_+1) = Dx_(i_)/Dx1(i)*Ex(i,Bi_(4)+1,k+1); 
	  }}
	  for(i_=0; i_<Nx-1; i_++){  
	   for(k_=0; k_<Nz; k_++){
	    i = ... ;
	    k = ... ;   		
	    Ez_(i_+1,0  ,k_) = Dz1_(k_)/Dz1(k)*Ez(i+1,Bi_(1)  ,k); 
            Ez_(i_+1,Ny_,k_) = Dz1_(k_)/Dz1(k)*Ez(i+1,Bi_(4)+1,k); 
	  }}

	  // z-faces
          for(i_=0; i_<Nx; i_++){  
	   for(j_=0; j_<Ny-1; j_++){
	    i = ... ;
	    j = ... ;               		
	    Ex_(i_,j_+1,0  ) = Dx_(i_)/Dx1(i)*Ex(i,j+1,Bi_(2)  ); 
            Ex_(i_,j_+1,Nz_) = Dx_(i_)/Dx1(i)*Ex(i,j+1,Bi_(5)+1); 
	  }}
	  for(i_=0; i_<Nx-1; i_++){  
	   for(j_=0; j_<Ny; j_++){ 
            i = ... ;
	    j = ... ;  		
	    Ey(i_+1,j_,0  ) = Dy1_(j_)/Dy1(j)*Ey(i+1,j,Bi_(2)  );
            Ey(i_+1,j_,Nz_) = Dy1_(j_)/Dy1(j)*Ey(i+1,j,Bi_(5)+1);
          }}
                                      
	}

}



void freeMemorySubgrids(){
	
}


/* COMMENTS:
(1) We use the conventional time-average discretization of the conduction current.
(2) E components are shifted by half a time step backward in time if they are solved implicitly. 
(3) Use the identity (I+A)^-1*(B-A) = -I + (I+A)^-1*(I+B) where A all off-diagonal entries are included in A.
(4) In order to preserve the structured memory of FDTD in the form of arrays, some extra flops are 
     done in the part of the main grid that overlaps with the subgrid. 
(5) The update coefficients of the coarse electric fields along the subgrid faces and edges are changed according
     to the finite-integration rule. The relative permeability of the coarse magnetic fields inside the overlapped
     part are set to infinity such that these fields remain zero-valued (PMC).   
(6) The cells next to the exterior PECs are not allowed to be subgridded because the time-difference auxiliary unknowns 
     are not defined in this PECs. This could be mittigated by implementing some if-statements. 
*/


