
#include "adhie.h"


void initAdhie1(Grid *g, PyObject *ipt){

	int s, i, j, k, l;    

	// Python input data
	PyObject *py_mg    = PyObject_GetAttrString( ipt, "maingrid" );
	PyObject *py_adhie = PyObject_GetAttrString( py_mg, "adhie" );
	PyObject *py_thr   = PyObject_GetAttrString( py_adhie, "thr" );
	PyObject *py_alpha = PyObject_GetAttrString( py_adhie, "alpha" );
	double dxt = PyFloat_AsDouble( PyList_GetItem( py_thr, 0 ) );       // step thresholds
	double dyt = PyFloat_AsDouble( PyList_GetItem( py_thr, 1 ) ); 	
	double dzt = PyFloat_AsDouble( PyList_GetItem( py_thr, 2 ) );           
	double alpha = 4.0*PyFloat_AsDouble(py_alpha);                      // adhie accuracy parameter 
	Py_DECREF(py_alpha);
	Py_DECREF(py_thr);
	Py_DECREF(py_adhie);
	Py_DECREF(py_mg);	

	// determine the number of slices or, equivalently, the number of tridiagonal blocks		
	Sx=(Dy1(0)<=dyt); for(j=1; j<Ny ; j++){ if(Dy1(j)<=dyt && Dy1(j-1)>dyt){ Sx++; }} 		
	Sy=(Dz1(0)<=dzt); for(k=1; k<Nz ; k++){ if(Dz1(k)<=dzt && Dz1(k-1)>dzt){ Sy++; }} 
	Sz=(Dx1(0)<=dxt); for(i=1; i<Nx ; i++){ if(Dx1(i)<=dxt && Dx1(i-1)>dxt){ Sz++; }} 

	// allocate memory for the slices
	my_alloc(slx, Sx, Slice);			
	my_alloc(sly, Sy, Slice);	
	my_alloc(slz, Sz, Slice);	

	// determine the start-index and size of each tridiagonal block 
	s=0; 	
	for(j=0; j<Ny; j++){ 
	  if(Dy1(j)<=dyt){ 	
	    Ix(s)=max(j-1,0); 
	    while(j+1<Ny && Dy1(j+1)<=dyt){ j++; } 
	    Nyh(s)=min(j+1-Ix(s),Ny-1-Ix(s));	    	 		     	
	    s++; 	                          }}

	s=0; 	
	for(k=0; k<Nz; k++){ 
	  if(Dz1(k)<=dzt){ 
	    Iy(s)=max(k-1,0); 	
	    while(k+1<Nz && Dz1(k+1)<=dzt){ k++; } 
	    Nzh(s)=min(k+1-Iy(s),Nz-1-Iy(s));	    		      	
            s++; 	                         }}

	s=0; 	
	for(i=0; i<Nx; i++){ 
	  if(Dx1(i)<=dxt){ 	
	    Iz(s)=max(i-1,0); 
	    while(i+1<Nx && Dx1(i+1)<=dxt){ i++; } 
	    Nxh(s)=min(i+1-Iz(s),Nx-1-Iz(s));	    		     	
            s++; 	                         }}

	for(s=0; s<Sx; s++){ Nye(s)=Nyh(s)+1-(Ix(s)==0 ? 1:0)-(Ix(s)+Nyh(s)==Ny-1 ? 1:0); }
	for(s=0; s<Sy; s++){ Nze(s)=Nzh(s)+1-(Iy(s)==0 ? 1:0)-(Iy(s)+Nzh(s)==Nz-1 ? 1:0); }
        for(s=0; s<Sz; s++){ Nxe(s)=Nxh(s)+1-(Iz(s)==0 ? 1:0)-(Iz(s)+Nxh(s)==Nx-1 ? 1:0); }

	// allocate memory for the update coefficients
	for(s=0; s<Sx; s++){ my_alloc(slx[s].qe, Nx*Nye(s)*(Nz-1)*3, double); my_alloc(slx[s].qh, (Nx-1)*Nyh(s)*Nz*3, double); }	
	for(s=0; s<Sy; s++){ my_alloc(sly[s].qe, (Nx-1)*Ny*Nze(s)*3, double); my_alloc(sly[s].qh, Nx*(Ny-1)*Nzh(s)*3, double); }	
	for(s=0; s<Sz; s++){ my_alloc(slz[s].qe, Nxe(s)*(Ny-1)*Nz*3, double); my_alloc(slz[s].qh, Nxh(s)*Ny*(Nz-1)*3, double); }

	// update coefficients for the Ex                                                
	for(s=0; s<Sx; s++){	
	  for(i=0; i<Nx; i++){  
	    for(k=0; k<Nz-1; k++){   
	      for(j=0; j<Nye(s); j++){
	        l=Ix(s)+j; 		      
	        Qex(s,i,j,k,0) = Cex(i,l,k)/alpha;     // subdiagonal
	        Qex(s,i,j,k,2) = Cex(i,l,k)/alpha;     // superdiagonal  	                  	
	}}}}

	// update coefficients for Ey
	for(s=0; s<Sy; s++){	
	  for(i=0; i<Nx-1; i++){  
	    for(j=0; j<Ny; j++){
	      for(k=0; k<Nze(s); k++){
		l=Iy(s)+k; 	
	        Qey(s,i,j,k,0) = Cey(i,j,l)/alpha;    
	        Qey(s,i,j,k,2) = Cey(i,j,l)/alpha;                    	
	}}}}
                                
	// update coefficients for Ez		
	for(s=0; s<Sz; s++){	 	
	  for(j=0; j<Ny-1; j++){
	    for(k=0; k<Nz; k++){
              for(i=0; i<Nxe(s); i++){  
	        l=Iz(s)+i;  
                Qez(s,i,j,k,0) = Cez(l,j,k)/alpha;    
	        Qez(s,i,j,k,2) = Cez(l,j,k)/alpha;                   	
	}}}}
	
	// update coefficients for Hx
	for(s=0; s<Sx; s++){ 	
	  for(i=0; i<Nx-1; i++){  
	    for(k=0; k<Nz; k++){                                     
	      for(j=0; j<Nyh(s); j++){   
		l = Ix(s)+j; 	        
		Qhx(s,i,j,k,0) = Chx(i,l,k)/alpha;      
	        Qhx(s,i,j,k,2) = Chx(i,l,k)/alpha;        
	}}}}

	// update coefficients for Hy
	for(s=0; s<Sy; s++){ 	
	  for(i=0; i<Nx; i++){  
	    for(j=0; j<Ny-1; j++){                                     
	      for(k=0; k<Nzh(s); k++){
		l=Iy(s)+k; 	        
		Qhy(s,i,j,k,0) = Chy(i,j,l)/alpha;      
	        Qhy(s,i,j,k,2) = Chy(i,j,l)/alpha;           
	}}}}

	// update coefficients for Hz
	for(s=0; s<Sz; s++){ 	
	  for(j=0; j<Ny; j++){   
	    for(k=0; k<Nz-1; k++){                                      
	      for(i=0; i<Nxh(s); i++){ 
		l=Iz(s)+i; 	        
		Qhz(s,i,j,k,0) = Chz(l,j,k)/alpha;      
	        Qhz(s,i,j,k,2) = Chz(l,j,k)/alpha;         
	}}}}

}



void initAdhie2(Grid *g){

	int s, i, j, k, l;  

	// update coefficients for the Ex                                                
	for(s=0; s<Sx; s++){	
	  for(i=0; i<Nx; i++){  
	    for(k=0; k<Nz-1; k++){   
	      j=0; l=Ix(s)+j;              
	      Qex(s,i,j,k,0) *= Chz(i,l  ,k);
	      Qex(s,i,j,k,2) *= Chz(i,l+1,k);
	      Qex(s,i,j,k,1) = 1.0-Qex(s,i,j,k,0)-Qex(s,i,j,k,2);   
	      Qex(s,i,j,k,1) = 1.0/Qex(s,i,j,k,1);   
	      Qex(s,i,j,k,2) *= Qex(s,i,j,k,1);	
	      for(j=1; j<Nye(s)-1; j++){
	        l=Ix(s)+j; 		      
	        Qex(s,i,j,k,0) *= Chz(i,l  ,k);                                           // subdiagonal (negative)
	        Qex(s,i,j,k,2) *= Chz(i,l+1,k);                                           // superdiagonal (negative)  
	        Qex(s,i,j,k,1) = 1.0-Qex(s,i,j,k,0)-Qex(s,i,j,k,2);                       // main diagonal (positive)
		Qex(s,i,j,k,1) = 1.0/(Qex(s,i,j,k,1)-Qex(s,i,j,k,0)*Qex(s,i,j-1,k,2));    // decomposition (Thomas algorithm) 
	        Qex(s,i,j,k,2) *= Qex(s,i,j,k,1);
	      }
	      j=Nye(s)-1; l=Ix(s)+j; 	
	      Qex(s,i,j,k,0) *= Chz(i,l  ,k);                    
	      Qex(s,i,j,k,2) *= Chz(i,l+1,k);                    
	      Qex(s,i,j,k,1) = 1.0-Qex(s,i,j,k,0)-Qex(s,i,j,k,2);                   		
              Qex(s,i,j,k,1) = 1.0/(Qex(s,i,j,k,1)-Qex(s,i,j,k,0)*Qex(s,i,j-1,k,2));   	                  	
	}}}

	// update coefficients for Ey
	for(s=0; s<Sy; s++){	
	  for(i=0; i<Nx-1; i++){  
	    for(j=0; j<Ny; j++){
	      k=0; l=Iy(s)+k;		      
	      Qey(s,i,j,k,0) *= Chx(i,j,l  );    
	      Qey(s,i,j,k,2) *= Chx(i,j,l+1);    
	      Qey(s,i,j,k,1) = 1.0-Qey(s,i,j,k,0)-Qey(s,i,j,k,2);
	      Qey(s,i,j,k,1) = 1.0/Qey(s,i,j,k,1); 
	      Qey(s,i,j,k,2) *= Qey(s,i,j,k,1);	
	      for(k=1; k<Nze(s)-1; k++){
		l=Iy(s)+k; 	
	        Qey(s,i,j,k,0) *= Chx(i,j,l  );    
	        Qey(s,i,j,k,2) *= Chx(i,j,l+1);    
	        Qey(s,i,j,k,1) = 1.0-Qey(s,i,j,k,0)-Qey(s,i,j,k,2);  
                Qey(s,i,j,k,1) = 1.0/(Qey(s,i,j,k,1)-Qey(s,i,j,k,0)*Qey(s,i,j,k-1,2));  
	        Qey(s,i,j,k,2) *= Qey(s,i,j,k,1); 
	      }
	      k=Nze(s)-1; l=Iy(s)+k; 
              Qey(s,i,j,k,0) *= Chx(i,j,l  );    
	      Qey(s,i,j,k,2) *= Chx(i,j,l+1);    
	      Qey(s,i,j,k,1) = 1.0-Qey(s,i,j,k,0)-Qey(s,i,j,k,2);   	
              Qey(s,i,j,k,1) = 1.0/(Qey(s,i,j,k,1)-Qey(s,i,j,k,0)*Qey(s,i,j,k-1,2));                  	
	}}}
                                
	// update coefficients for Ez		
	for(s=0; s<Sz; s++){	 	
	  for(j=0; j<Ny-1; j++){
	    for(k=0; k<Nz; k++){
	      i=0; l=Iz(s)+i;  
              Qez(s,i,j,k,0) *= Chy(l  ,j,k);    
	      Qez(s,i,j,k,2) *= Chy(l+1,j,k);   
	      Qez(s,i,j,k,1) = 1.0-Qez(s,i,j,k,0)-Qez(s,i,j,k,2);
	      Qez(s,i,j,k,1) = 1.0/Qez(s,i,j,k,1);  	 	              
	      Qez(s,i,j,k,2) *= Qez(s,i,j,k,1);	
              for(i=1; i<Nxe(s)-1; i++){  
	        l=Iz(s)+i;  
                Qez(s,i,j,k,0) *= Chy(l  ,j,k);    
	        Qez(s,i,j,k,2) *= Chy(l+1,j,k);   
	        Qez(s,i,j,k,1) = 1.0-Qez(s,i,j,k,0)-Qez(s,i,j,k,2);  	
                Qez(s,i,j,k,1) = 1.0/(Qez(s,i,j,k,1)-Qez(s,i,j,k,0)*Qez(s,i-1,j,k,2)); 
	        Qez(s,i,j,k,2) *= Qez(s,i,j,k,1);
	      }
              i=Nxe(s)-1; l=Iz(s)+i;  
              Qez(s,i,j,k,0) *= Chy(l  ,j,k);    
	      Qez(s,i,j,k,2) *= Chy(l+1,j,k);   
	      Qez(s,i,j,k,1) = 1.0-Qez(s,i,j,k,0)-Qez(s,i,j,k,2); 	
              Qez(s,i,j,k,1) = 1.0/(Qez(s,i,j,k,1)-Qez(s,i,j,k,0)*Qez(s,i-1,j,k,2));                	
	}}}

	// update coefficients for Hx
	for(s=0; s<Sx; s++){ 	
	  for(i=0; i<Nx-1; i++){  
	    for(k=0; k<Nz; k++){       
	      j=0; l=Ix(s)+j;
	      Qhx(s,i,j,k,0) *= (l==0 ? 0.0 : 2.0*Cez(i,l-1,k)/(2.0+Sez(i,l-1,k)));  
	      Qhx(s,i,j,k,2) *=               2.0*Cez(i,l  ,k)/(2.0+Sez(i,l  ,k));       
	      Qhx(s,i,j,k,1) = 1.0-Qhx(s,i,j,k,0)-Qhx(s,i,j,k,2);
	      Qhx(s,i,j,k,1) = 1.0/Qhx(s,i,j,k,1); 	
	      Qhx(s,i,j,k,2) *= Qhx(s,i,j,k,1);                                 
	      for(j=1; j<Nyh(s)-1; j++){   
		l = Ix(s)+j; 	        
		Qhx(s,i,j,k,0) *= 2.0*Cez(i,l-1,k)/(2.0+Sez(i,l-1,k));      
	        Qhx(s,i,j,k,2) *= 2.0*Cez(i,l  ,k)/(2.0+Sez(i,l  ,k));        
	        Qhx(s,i,j,k,1) = 1.0-Qhx(s,i,j,k,0)-Qhx(s,i,j,k,2);                             
	        Qhx(s,i,j,k,1) = 1.0/(Qhx(s,i,j,k,1)-Qhx(s,i,j,k,0)*Qhx(s,i,j-1,k,2));              
                Qhx(s,i,j,k,2) *= Qhx(s,i,j,k,1);                    
	        }
	      j=Nyh(s)-1; l=Ix(s)+j; 	
	      Qhx(s,i,j,k,0) *=                  2.0*Cez(i,l-1,k)/(2.0+Sez(i,l-1,k));      
	      Qhx(s,i,j,k,2) *= (l==Ny-1 ? 0.0 : 2.0*Cez(i,l  ,k)/(2.0+Sez(i,l  ,k)));
	      Qhx(s,i,j,k,1) = 1.0-Qhx(s,i,j,k,0)-Qhx(s,i,j,k,2);
              Qhx(s,i,j,k,1) = 1.0/(Qhx(s,i,j,k,1)-Qhx(s,i,j,k,0)*(j==0 ? 0.0 : Qhx(s,i,j-1,k,2))); 		          
	}}}

	// update coefficients for Hy
	for(s=0; s<Sy; s++){ 	
	  for(i=0; i<Nx; i++){  
	    for(j=0; j<Ny-1; j++){           
	      k=0; l=Iy(s)+k;
	      Qhy(s,i,j,k,0) *= (l==0 ? 0.0 : 2.0*Cex(i,j,l-1)/(2.0+Sex(i,j,l-1)));  
	      Qhy(s,i,j,k,2) *=               2.0*Cex(i,j,l  )/(2.0+Sex(i,j,l  ));      
	      Qhy(s,i,j,k,1) = 1.0-Qhy(s,i,j,k,0)-Qhy(s,i,j,k,2);
	      Qhy(s,i,j,k,1) = 1.0/Qhy(s,i,j,k,1);      	
	      Qhy(s,i,j,k,2) *= Qhy(s,i,j,k,1);                             
	      for(k=1; k<Nzh(s)-1; k++){
		l=Iy(s)+k; 	        
		Qhy(s,i,j,k,0) *= 2.0*Cex(i,j,l-1)/(2.0+Sex(i,j,l-1));      
	        Qhy(s,i,j,k,2) *= 2.0*Cex(i,j,l  )/(2.0+Sex(i,j,l  ));         
	        Qhy(s,i,j,k,1) = 1.0-Qhy(s,i,j,k,0)-Qhy(s,i,j,k,2);
		Qhy(s,i,j,k,1) = 1.0/(Qhy(s,i,j,k,1)-Qhy(s,i,j,k,0)*Qhy(s,i,j,k-1,2)); 
                Qhy(s,i,j,k,2) *= Qhy(s,i,j,k,1);                              
	        }
	      k=Nzh(s)-1; l=Iy(s)+k; 	
	      Qhy(s,i,j,k,0) *=                  2.0*Cex(i,j,l-1)/(2.0+Sex(i,j,l-1));      
	      Qhy(s,i,j,k,2) *= (l==Nz-1 ? 0.0 : 2.0*Cex(i,j,l  )/(2.0+Sex(i,j,l  )));         
	      Qhy(s,i,j,k,1) = 1.0-Qhy(s,i,j,k,0)-Qhy(s,i,j,k,2);
	      Qhy(s,i,j,k,1) = 1.0/(Qhy(s,i,j,k,1)-Qhy(s,i,j,k,0)*(k==0 ? 0.0 : Qhy(s,i,j,k-1,2)));	  
	}}}

	// update coefficients for Hz
	for(s=0; s<Sz; s++){ 	
	  for(j=0; j<Ny; j++){   
	    for(k=0; k<Nz-1; k++){     
	      i=0; l=Iz(s)+i;
	      Qhz(s,i,j,k,0) *= (l==0 ? 0.0 : 2.0*Cey(l-1,j,k)/(2.0+Sey(l-1,j,k)));
	      Qhz(s,i,j,k,2) *=               2.0*Cey(l  ,j,k)/(2.0+Sey(l  ,j,k));      
              Qhz(s,i,j,k,1) = 1.0-Qhz(s,i,j,k,0)-Qhz(s,i,j,k,2); 
	      Qhz(s,i,j,k,1) = 1.0/Qhz(s,i,j,k,1);  	
	      Qhz(s,i,j,k,2) *= Qhz(s,i,j,k,1);                                  
	      for(i=1; i<Nxh(s)-1; i++){ 
		l=Iz(s)+i; 	        
		Qhz(s,i,j,k,0) *= 2.0*Cey(l-1,j,k)/(2.0+Sey(l-1,j,k));      
	        Qhz(s,i,j,k,2) *= 2.0*Cey(l  ,j,k)/(2.0+Sey(l  ,j,k));         
	        Qhz(s,i,j,k,1) = 1.0-Qhz(s,i,j,k,0)-Qhz(s,i,j,k,2); 
                Qhz(s,i,j,k,1) = 1.0/(Qhz(s,i,j,k,1)-Qhz(s,i,j,k,0)*Qhz(s,i-1,j,k,2)); 
                Qhz(s,i,j,k,2) *= Qhz(s,i,j,k,1);                             
	        }
	      i=Nxh(s)-1; l=Iz(s)+i; 	
	      Qhz(s,i,j,k,0) *=                  2.0*Cey(l-1,j,k)/(2.0+Sey(l-1,j,k));    
              Qhz(s,i,j,k,2) *= (l==Nx-1 ? 0.0 : 2.0*Cey(l  ,j,k)/(2.0+Sey(l  ,j,k))); 		     
	      Qhz(s,i,j,k,1) = 1.0-Qhz(s,i,j,k,0)-Qhz(s,i,j,k,2); 
	      Qhz(s,i,j,k,1) = 1.0/(Qhz(s,i,j,k,1)-Qhz(s,i,j,k,0)*(i==0 ? 0.0 : Qhz(s,i-1,j,k,2))); 
	}}}

}




void updateDiffEadhie(Grid *g){

	int s, i, j, k, l;
	
	for(s=0; s<Sx; s++){	
	  for (i=0; i<Nx; i++){  
	    for (k=0; k<Nz-1; k++){
	      DiffEx(i,Ix(s),k) *= Qex(s,i,0,k,1);
	      for (j=1       ; j<Nye(s); j++){ l=Ix(s)+j; DiffEx(i,l,k) = (DiffEx(i,l,k)-Qex(s,i,j,k,0)*DiffEx(i,l-1,k))*Qex(s,i,j,k,1); }   // forward sweep                     
	      for (j=Nye(s)-2; j>=0    ; j--){ l=Ix(s)+j; DiffEx(i,l,k) -= Qex(s,i,j,k,2)*DiffEx(i,l+1,k); }                                 // backward sweep 	
	}}}

	for(s=0; s<Sy; s++){	
	  for (i=0; i<Nx-1; i++){  
	    for (j=0; j<Ny; j++){
	      DiffEy(i,j,Iy(s)) *= Qey(s,i,j,0,1);
	      for (k=1       ; k<Nze(s); k++){ l=Iy(s)+k; DiffEy(i,j,l) = (DiffEy(i,j,l)-Qey(s,i,j,k,0)*DiffEy(i,j,l-1))*Qey(s,i,j,k,1); }                       
	      for (k=Nze(s)-2; k>=0    ; k--){ l=Iy(s)+k; DiffEy(i,j,l) -= Qey(s,i,j,k,2)*DiffEy(i,j,l+1); }                                 	
	}}}

	for(s=0; s<Sz; s++){	
	  for (j=0; j<Ny-1; j++){  
	    for (k=0; k<Nz; k++){
	      DiffEz(Iz(s),j,k) *= Qez(s,0,j,k,1);
	      for (i=1       ; i<Nxe(s); i++){ l=Iz(s)+i; DiffEz(l,j,k) = (DiffEz(l,j,k)-Qez(s,i,j,k,0)*DiffEz(l-1,j,k))*Qez(s,i,j,k,1); }                       
	      for (i=Nxe(s)-2; i>=0    ; i--){ l=Iz(s)+i; DiffEz(l,j,k) -= Qez(s,i,j,k,2)*DiffEz(l+1,j,k); }                                 	
	}}}

}


void updateDiffHadhie(Grid *g){

	int s, i, j, k, l;

	for(s=0; s<Sx; s++){	
	  for (i=0; i<Nx-1; i++){  
	    for (k=0; k<Nz; k++){
	      DiffHx(i,Ix(s),k) *= Qhx(s,i,0,k,1);	
              for (j=1       ; j<Nyh(s); j++){ l=Ix(s)+j; DiffHx(i,l,k) = (DiffHx(i,l,k)-Qhx(s,i,j,k,0)*DiffHx(i,l-1,k))*Qhx(s,i,j,k,1); }   // forward sweep                         
	      for (j=Nyh(s)-2; j>=0    ; j--){ l=Ix(s)+j; DiffHx(i,l,k) -= Qhx(s,i,j,k,2)*DiffHx(i,l+1,k); }                                 // backward sweep 	          
	}}}

	for(s=0; s<Sy; s++){	
	  for (i=0; i<Nx; i++){  
	    for (j=0; j<Ny-1; j++){
	      DiffHy(i,j,Iy(s)) *= Qhy(s,i,j,0,1);	
              for (k=1       ; k<Nzh(s); k++){ l=Iy(s)+k; DiffHy(i,j,l) = (DiffHy(i,j,l)-Qhy(s,i,j,k,0)*DiffHy(i,j,l-1))*Qhy(s,i,j,k,1); }                           
	      for (k=Nzh(s)-2; k>=0    ; k--){ l=Iy(s)+k; DiffHy(i,j,l) -= Qhy(s,i,j,k,2)*DiffHy(i,j,l+1); }                                 	          
	}}}

	for(s=0; s<Sz; s++){	
	  for (j=0; j<Ny; j++){  
	    for (k=0; k<Nz-1; k++){
	      DiffHz(Iz(s),j,k) *= Qhz(s,0,j,k,1);	
              for (i=1       ; i<Nxh(s); i++){ l=Iz(s)+i; DiffHz(l,j,k) = (DiffHz(l,j,k)-Qhz(s,i,j,k,0)*DiffHz(l-1,j,k))*Qhz(s,i,j,k,1); }                           
	      for (i=Nxh(s)-2; i>=0    ; i--){ l=Iz(s)+i; DiffHz(l,j,k) -= Qhz(s,i,j,k,2)*DiffHz(l+1,j,k); }                                 	          
	}}}

}



void freeMemoryAdhie(){
	int s; 	
	for (s=Sz-1; s>=0; s--){ free(slz[s].qh); free(slz[s].qe); } 
	for (s=Sy-1; s>=0; s--){ free(sly[s].qh); free(sly[s].qe); }	
	for (s=Sx-1; s>=0; s--){ free(slx[s].qh); free(slx[s].qe); }  
	free(slz); free(sly); free(slx);   
}








/* COMMENTS:
(1) ADI requires a total of 3 additional coefficients per update  
(2) Avoid divisions in the time stepping loop (use sums and multiplications as much as possible).
(3) The tridiagonal solver can be parallellized (see Numerical recipies in C p.57).
(4) By construction, the tridiagonal matrices are diagonally dominant: |main diag|>|subdiag|+|superdiag|, 
      such that the Thomas algorithm cannot encounter a zero pivot.
(5) For lossy media, the non-trivial self-term is computed using the matrix identity:
      (A+B)^-1 * (A-B) = I - 2 * (A+B)^-1 * B .      
    Left-multiply both sides by (A+B) to verify this identity. 
    Hence,
      E = (A+B)^-1 * [(A-B)*E + DiffE + J]  
    reduces to 
      E = E + (A+B)^-1 * [DiffE + J - 2*B*E] 
(6) Only for implicitization in the z-dimension, the loops are nested in the most efficient order.  
      Implicitization in the x-dimension should be avoided.  	   	
(7) The real-stretch PML-parameter kappa only occurs once in the ADI perturbation term instead of twice (doi:10.1109/TEMC.2012.2198067). 
      This is resolved by splitting the ADHIE initialization into two subroutines. One is executed before and the other 
      after the pml initialization. 
(8) The output of the tridiagonal solver has been compared to Lapack's dgtsv with success.   
*/






