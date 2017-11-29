
#include "pml.h"


void initPml(Grid *g, PyObject *ipt){	

	// sigma2t number of PML layers  
	PyObject *py_mg   = PyObject_GetAttrString( ipt, "maingrid" );	
	PyObject *py_Npml = PyObject_GetAttrString( py_mg, "Npml" ); 
	for (int i=0;i<6;i++){ Npml(i) = PyInt_AS_LONG(PyList_GET_ITEM( py_Npml, (Py_ssize_t) i)); }	
	Py_DECREF(py_Npml);

	// dynamically allocate memory for the PMLs
	Ntan(0,0)=Ny; Ntan(0,1)=Nz; my_alloc_pml(0); 
	Ntan(1,0)=Ny; Ntan(1,1)=Nz; my_alloc_pml(1);
	Ntan(2,0)=Nz; Ntan(2,1)=Nx; my_alloc_pml(2); 
	Ntan(3,0)=Nz; Ntan(3,1)=Nx; my_alloc_pml(3);
	Ntan(4,0)=Nx; Ntan(4,1)=Ny; my_alloc_pml(4);
	Ntan(5,0)=Nx; Ntan(5,1)=Ny; my_alloc_pml(5);
	
	// Python objects: maximum real and CFS stretch 	
	PyObject *py_Kmax = PyObject_GetAttrString( py_mg, "Kmax" ); 
	PyObject *py_Amax = PyObject_GetAttrString( py_mg, "Amax" );

	// primary and dual stretch quantities 	  
	double polygrad1, polygrad2, kappa1, kappa2, sigma1, sigma2, alpha1, alpha2;   

	// three spatial indices, PML index, translated i index
	int i, j, k, p, it; 
	                              
	p=0; // left x-pml (i=x, i=y, k=z)		
	for (i=0; i<Npml(p); i++){                        	

	    // polynomial grading  
	    polygrad1 = pow((double)(Npml(p)-i    )/Npml(p),m_pg); 
	    polygrad2 = pow((double)(Npml(p)-i-0.5)/Npml(p),m_pg);         
	    
            // real-stretch parameter
            kappa1 = 1.0+(Kmax(p)-1.0)*polygrad1; 
            kappa2 = 1.0+(Kmax(p)-1.0)*polygrad2;                         
	    
            // imaginary-stretch parameter (DOI:10.1515/jee-2017-0006)
            sigma1 = 0.8*(m_pg+1.0)/Dx1(i)*polygrad1;  
            sigma2 = 0.8*(m_pg+1.0)/Dx2(i)*polygrad2;                     
	    
            // linearly graded complex-frequency shift 
            alpha1 = Z0*Amax(p)*(double)(i+1.0)/Npml(p)+sigma1/kappa1; 
            alpha2 = Z0*Amax(p)*(double)(i+0.5)/Npml(p)+sigma2/kappa2;   

            // first coefficient of the auxiliary equation (convolutional formulation) 
            Ae(p,i,0) = exp(-alpha2*Dtau);                                 
       	    Ah(p,i,0) = exp(-alpha1*Dtau);                                 

            // second coefficient of the auxiliary equation (convolutional formulation)
            Ae(p,i,1) = sigma2/(kappa2*alpha2)*(Ae(p,i,0)-1.0);         
       	    Ah(p,i,1) = sigma1/(kappa1*alpha1)*(Ah(p,i,0)-1.0);          		

            // change standard FDTD coefficient in front of the curl  
            for (j=0; j<Ntan(p,0); j++){                                 
 	      for (k=0; k<Ntan(p,1)-1; k++){                 	                                     
                Cey(i,j,k) /= kappa2;		                     
                Chz(i,j,k) /= kappa1;				      
	    }}    
            for (j=0; j<Ntan(p,0)-1; j++){
	      for (k=0; k<Ntan(p,1); k++){ 
	        Cez(i,j,k) /= kappa2;	
		Chy(i,j,k) /= kappa1;			      
	    }}
            for (j=0; j<Ntan(p,0)-1; j++){                                 
 	      for (k=0; k<Ntan(p,1)-1; k++){                 	                                     
                Cex(i,j,k) *= kappa1;		                     			      
	    }}
            for (j=0; j<Ntan(p,0); j++){                                 
 	      for (k=0; k<Ntan(p,1); k++){                 	                                     
                Chx(i,j,k) *= kappa2;		                     				      
	    }}
  
        }
	   
	p=1; // right x-pml (i=x, i=y, k=z)
	for (i=0; i<Npml(p); i++){                       
	    
	    it = Nx-Npml(p)+i; 

	    polygrad1 = pow((double)(i+1.0)/Npml(p),m_pg); 
	    polygrad2 = pow((double)(i+0.5)/Npml(p),m_pg);           
	    
	    kappa1 = 1.0+(Kmax(p)-1.0)*polygrad1; 
            kappa2 = 1.0+(Kmax(p)-1.0)*polygrad2;                        

	    sigma1 = 0.8*(m_pg+1.0)/Dx1(it  )*polygrad1;
            sigma2 = 0.8*(m_pg+1.0)/Dx2(it-1)*polygrad2;         
	    
	    alpha1 = Z0*Amax(p)*(double)(Npml(p)-i    )/Npml(p)+sigma1/kappa1;
            alpha2 = Z0*Amax(p)*(double)(Npml(p)-i-0.5)/Npml(p)+sigma2/kappa2;    
            
            Ae(p,i,0) = exp(-alpha2*Dtau);                        
       	    Ah(p,i,0) = exp(-alpha1*Dtau); 
	    
            Ae(p,i,1) = sigma2/(kappa2*alpha2)*(Ae(p,i,0)-1.0);        
       	    Ah(p,i,1) = sigma1/(kappa1*alpha1)*(Ah(p,i,0)-1.0); 			
            
            for (j=0; j<Ntan(p,0); j++){                     
 	      for (k=0; k<Ntan(p,1)-1; k++){                 	 	        
		Cey(it-1,j,k) /= kappa2;	            
	        Chz(it  ,j,k) /= kappa1;				      
	    }}
            for (j=0; j<Ntan(p,0)-1; j++){
	      for (k=0; k<Ntan(p,1); k++){
	        Cez(it-1,j,k) /= kappa2;	
	        Chy(it  ,j,k) /= kappa1;			      
	    }}
            for (j=0; j<Ntan(p,0)-1; j++){                                 
 	      for (k=0; k<Ntan(p,1)-1; k++){                 	                                     
                Cex(it  ,j,k) *= kappa1;		                     			      
	    }}
            for (j=0; j<Ntan(p,0); j++){                                 
 	      for (k=0; k<Ntan(p,1); k++){                 	                                     
                Chx(it-1,j,k) *= kappa2;		                     				      
	    }}           

        }	
                         
	p=2; // left y-pml (i=y, j=z, k=x) 		
	for (i=0; i<Npml(p); i++){                        	

	    polygrad1 = pow((double)(Npml(p)-i    )/Npml(p),m_pg); 
	    polygrad2 = pow((double)(Npml(p)-i-0.5)/Npml(p),m_pg);        
	    
            kappa1 = 1.0+(Kmax(p)-1.0)*polygrad1; 
            kappa2 = 1.0+(Kmax(p)-1.0)*polygrad2;                         
	    
	    sigma1 = 0.8*(m_pg+1.0)/Dy1(i)*polygrad1; 
            sigma2 = 0.8*(m_pg+1.0)/Dy2(i)*polygrad2;            
	    
	    alpha1 = Z0*Amax(p)*(double)(i+1.0)/Npml(p)+sigma1/kappa1; 
            alpha2 = Z0*Amax(p)*(double)(i+0.5)/Npml(p)+sigma2/kappa2;          

            Ae(p,i,0) = exp(-alpha2*Dtau);                                 
       	    Ah(p,i,0) = exp(-alpha1*Dtau); 

            Ae(p,i,1) = sigma2/(kappa2*alpha2)*(Ae(p,i,0)-1.0);  
       	    Ah(p,i,1) = sigma1/(kappa1*alpha1)*(Ah(p,i,0)-1.0);			
	    	
            for (j=0; j<Ntan(p,0); j++){                                 
 	      for (k=0; k<Ntan(p,1)-1; k++){                 	                               
                Cez(k,i,j) /= kappa2;		        
                Chx(k,i,j) /= kappa1;				      
	    }}    
            for (j=0; j<Ntan(p,0)-1; j++){
	      for (k=0; k<Ntan(p,1); k++){
	        Cex(k,i,j) /= kappa2;	
		Chz(k,i,j) /= kappa1;			      
	    }} 
	    for (j=0; j<Ntan(p,0)-1; j++){                                 
 	      for (k=0; k<Ntan(p,1)-1; k++){                 	                               
                Cey(k,i,j) *= kappa1;		        		      
	    }} 
            for (j=0; j<Ntan(p,0); j++){                                 
 	      for (k=0; k<Ntan(p,1); k++){                 	                               
                Chy(k,i,j) *= kappa2;		        		      
	    }} 
 
        }
	                                                       
	p=3; // right y-pml (i=y, j=z, k=x)  
	for (i=0; i<Npml(p); i++){                        

	    it = Ny-Npml(p)+i;        

	    polygrad1 = pow((double)(i+1.0)/Npml(p),m_pg); 
	    polygrad2 = pow((double)(i+0.5)/Npml(p),m_pg);           
	    
	    kappa1 = 1.0+(Kmax(p)-1.0)*polygrad1; 
            kappa2 = 1.0+(Kmax(p)-1.0)*polygrad2;                        
	    
            sigma1 = 0.8*(m_pg+1.0)/Dy1(it  )*polygrad1;
            sigma2 = 0.8*(m_pg+1.0)/Dy2(it-1)*polygrad2;
	    
            alpha1 = Z0*Amax(p)*(double)(Npml(p)-i    )/Npml(p)+sigma1/kappa1;
            alpha2 = Z0*Amax(p)*(double)(Npml(p)-i-0.5)/Npml(p)+sigma2/kappa2;    
            
            Ae(p,i,0) = exp(-alpha2*Dtau);                        
       	    Ah(p,i,0) = exp(-alpha1*Dtau); 
	    
            Ae(p,i,1) = sigma2/(kappa2*alpha2)*(Ae(p,i,0)-1.0);        
       	    Ah(p,i,1) = sigma1/(kappa1*alpha1)*(Ah(p,i,0)-1.0); 			
            
            for (j=0; j<Ntan(p,0); j++){             
 	      for (k=0; k<Ntan(p,1)-1; k++){                 	 	        
		Cez(k,it-1,j) /= kappa2;	            
	        Chx(k,it  ,j) /= kappa1;				      
	    }}
            for (j=0; j<Ntan(p,0)-1; j++){
	      for (k=0; k<Ntan(p,1); k++){
	        Cex(k,it-1,j) /= kappa2;	
	        Chz(k,it  ,j) /= kappa1;			      
	    }}
	    for (j=0; j<Ntan(p,0)-1; j++){                                 
 	      for (k=0; k<Ntan(p,1)-1; k++){                 	                               
                Cey(k,it  ,j) *= kappa1;		        		      
	    }} 
            for (j=0; j<Ntan(p,0); j++){                                 
 	      for (k=0; k<Ntan(p,1); k++){                 	                               
                Chy(k,it-1,j) *= kappa2;		        		      
	    }} 

        }	
                            
	p=4; // left z-pml (i=z, j=x, k=y)  	
	for (i=0; i<Npml(p); i++){                        	

            polygrad1 = pow((double)(Npml(p)-i    )/Npml(p),m_pg); 
	    polygrad2 = pow((double)(Npml(p)-i-0.5)/Npml(p),m_pg);        
	    
            kappa1 = 1.0+(Kmax(p)-1.0)*polygrad1; 
            kappa2 = 1.0+(Kmax(p)-1.0)*polygrad2;                         
	    
	    sigma1 = 0.8*(m_pg+1.0)/Dz1(i)*polygrad1; 
            sigma2 = 0.8*(m_pg+1.0)/Dz2(i)*polygrad2;            
	    
	    alpha1 = Z0*Amax(p)*(double)(i+1.0)/Npml(p)+sigma1/kappa1; 
            alpha2 = Z0*Amax(p)*(double)(i+0.5)/Npml(p)+sigma2/kappa2;          

            Ae(p,i,0) = exp(-alpha2*Dtau);                                 
       	    Ah(p,i,0) = exp(-alpha1*Dtau); 

            Ae(p,i,1) = sigma2/(kappa2*alpha2)*(Ae(p,i,0)-1.0);  
       	    Ah(p,i,1) = sigma1/(kappa1*alpha1)*(Ah(p,i,0)-1.0);			
	    	
            for (j=0; j<Ntan(p,0); j++){                                 
 	      for (k=0; k<Ntan(p,1)-1; k++){                 	                                
                Cex(j,k,i) /= kappa2;		        
                Chy(j,k,i) /= kappa1;				      
	    }}    
            for (j=0; j<Ntan(p,0)-1; j++){
	      for (k=0; k<Ntan(p,1); k++){
	        Cey(j,k,i) /= kappa2;	
		Chx(j,k,i) /= kappa1;			      
	    }}  
	    for (j=0; j<Ntan(p,0)-1; j++){                                 
 	      for (k=0; k<Ntan(p,1)-1; k++){                 	                                
                Cez(j,k,i) *= kappa1;		        	      
	    }} 
	    for (j=0; j<Ntan(p,0); j++){                                 
 	      for (k=0; k<Ntan(p,1); k++){                 	                                
                Chz(j,k,i) *= kappa2;		        	      
	    }} 

        }
	                                              
	p=5; // right z-pml (i=z, j=x, k=y)           
	for (i=0; i<Npml(p); i++){                       
	    
	    it = Nz-Npml(p)+i; 

	    polygrad1 = pow((double)(i+1.0)/Npml(p),m_pg); 
	    polygrad2 = pow((double)(i+0.5)/Npml(p),m_pg);           
	    
            kappa1 = 1.0+(Kmax(p)-1.0)*polygrad1; 
            kappa2 = 1.0+(Kmax(p)-1.0)*polygrad2;                        
	    
	    sigma1 = 0.8*(m_pg+1.0)/Dz1(it  )*polygrad1;
            sigma2 = 0.8*(m_pg+1.0)/Dz2(it-1)*polygrad2;         
	    
	    alpha1 = Z0*Amax(p)*(double)(Npml(p)-i    )/Npml(p)+sigma1/kappa1;
            alpha2 = Z0*Amax(p)*(double)(Npml(p)-i-0.5)/Npml(p)+sigma2/kappa2;    
            
            Ae(p,i,0) = exp(-alpha2*Dtau);                        
       	    Ah(p,i,0) = exp(-alpha1*Dtau); 
	    
            Ae(p,i,1) = sigma2/(kappa2*alpha2)*(Ae(p,i,0)-1.0);        
       	    Ah(p,i,1) = sigma1/(kappa1*alpha1)*(Ah(p,i,0)-1.0); 			
            
            for (j=0; j<Ntan(p,0); j++){                     
 	      for (k=0; k<Ntan(p,1)-1; k++){                 	 	        
		Cex(j,k,it-1) /= kappa2;	            
	        Chy(j,k,it  ) /= kappa1;				      
	    }}
            for (j=0; j<Ntan(p,0)-1; j++){
	      for (k=0; k<Ntan(p,1); k++){
	        Cey(j,k,it-1) /= kappa2;	
	        Chx(j,k,it  ) /= kappa1;			      
	    }}
	    for (j=0; j<Ntan(p,0)-1; j++){                     
 	      for (k=0; k<Ntan(p,1)-1; k++){                 	 	        
		Cez(j,k,it  ) *= kappa1;	            				      
	    }}
	    for (j=0; j<Ntan(p,0); j++){                     
 	      for (k=0; k<Ntan(p,1); k++){                 	 	        
		Chz(j,k,it-1) *= kappa2;	            				      
	    }}

        }	
	

	// reference counting
	Py_DECREF(py_Kmax); 	
	Py_DECREF(py_Amax);	
	Py_DECREF(py_mg);

}



void updateDiffEpml(Grid *g){   
	
	int i, j, k, p, it;  // three spatial indices, PML index, translated i index

	p=0; // left x-pml	
	for (i=0; i<Npml(p); i++){  
	  for (j=0; j<Ntan(p,0); j++){
	    for (k=0; k<Ntan(p,1)-1; k++){ 		
	      E1(p,i,j,k)    = Ae(p,i,0)*E1(p,i,j,k) + Ae(p,i,1)*(Hz(i+1,j,k)-Hz(i,j,k));
	      DiffEy(i,j,k) -= Cey(i,j,k)*E1(p,i,j,k);
	  }}
	  for (j=0; j<Ntan(p,0)-1; j++){
	    for (k=0; k<Ntan(p,1); k++){ 		
	      E2(p,i,j,k)    = Ae(p,i,0)*E2(p,i,j,k) + Ae(p,i,1)*(Hy(i+1,j,k)-Hy(i,j,k));
	      DiffEz(i,j,k) += Cez(i,j,k)*E2(p,i,j,k);
	  }}
	}

	p=1; // right x-pml	
	for (i=0; i<Npml(p); i++){  
	  it = Nx-1-Npml(p)+i;
	  for (j=0; j<Ntan(p,0); j++){
	    for (k=0; k<Ntan(p,1)-1; k++){ 		
	      E1(p,i,j,k)     = Ae(p,i,0)*E1(p,i,j,k) + Ae(p,i,1)*(Hz(it+1,j,k)-Hz(it,j,k));
	      DiffEy(it,j,k) -= Cey(it,j,k)*E1(p,i,j,k);
	  }}
	  for (j=0; j<Ntan(p,0)-1; j++){
	    for (k=0; k<Ntan(p,1); k++){ 		
	      E2(p,i,j,k)     = Ae(p,i,0)*E2(p,i,j,k) + Ae(p,i,1)*(Hy(it+1,j,k)-Hy(it,j,k));
	      DiffEz(it,j,k) += Cez(it,j,k)*E2(p,i,j,k);
	  }}
	}

	p=2; // left y-pml	
	for (i=0; i<Npml(p); i++){  
	  for (j=0; j<Ntan(p,0); j++){                                                          
	    for (k=0; k<Ntan(p,1)-1; k++){ 		
	      E1(p,i,j,k)    = Ae(p,i,0)*E1(p,i,j,k) + Ae(p,i,1)*(Hx(k,i+1,j)-Hx(k,i,j));
	      DiffEz(k,i,j) -= Cez(k,i,j)*E1(p,i,j,k);
	  }}
	  for (j=0; j<Ntan(p,0)-1; j++){
	    for (k=0; k<Ntan(p,1); k++){ 		
	      E2(p,i,j,k)    = Ae(p,i,0)*E2(p,i,j,k) + Ae(p,i,1)*(Hz(k,i+1,j)-Hz(k,i,j));
	      DiffEx(k,i,j) += Cex(k,i,j)*E2(p,i,j,k);
	  }}
	}

	p=3; // right y-pml
	for (i=0; i<Npml(p); i++){  
	  it = Ny-1-Npml(p)+i; 
	  for (j=0; j<Ntan(p,0); j++){ 
	    for (k=0; k<Ntan(p,1)-1; k++){ 		
	      E1(p,i,j,k)     = Ae(p,i,0)*E1(p,i,j,k) + Ae(p,i,1)*(Hx(k,it+1,j)-Hx(k,it,j));
	      DiffEz(k,it,j) -= Cez(k,it,j)*E1(p,i,j,k);
	  }}
	  for (j=0; j<Ntan(p,0)-1; j++){
	    for (k=0; k<Ntan(p,1); k++){ 		
	      E2(p,i,j,k)     = Ae(p,i,0)*E2(p,i,j,k) + Ae(p,i,1)*(Hz(k,it+1,j)-Hz(k,it,j));
	      DiffEx(k,it,j) += Cex(k,it,j)*E2(p,i,j,k);
	  }}
	}

	p=4; // left z-pml	
	for (i=0; i<Npml(p); i++){  
	  for (j=0; j<Ntan(p,0); j++){
	    for (k=0; k<Ntan(p,1)-1; k++){ 		
	      E1(p,i,j,k)    = Ae(p,i,0)*E1(p,i,j,k) + Ae(p,i,1)*(Hy(j,k,i+1)-Hy(j,k,i));
	      DiffEx(j,k,i) -= Cex(j,k,i)*E1(p,i,j,k);
	  }}
	  for (j=0; j<Ntan(p,0)-1; j++){
	    for (k=0; k<Ntan(p,1); k++){ 		
	      E2(p,i,j,k)    = Ae(p,i,0)*E2(p,i,j,k) + Ae(p,i,1)*(Hx(j,k,i+1)-Hx(j,k,i));
	      DiffEy(j,k,i) += Cey(j,k,i)*E2(p,i,j,k);
	  }}
	}

	p=5; // right z-pml
	for (i=0; i<Npml(p); i++){  
	  it = Nz-1-Npml(p)+i;
	  for (j=0; j<Ntan(p,0); j++){
	    for (k=0; k<Ntan(p,1)-1; k++){ 		
	      E1(p,i,j,k)     = Ae(p,i,0)*E1(p,i,j,k) + Ae(p,i,1)*(Hy(j,k,it+1)-Hy(j,k,it));
	      DiffEx(j,k,it) -= Cex(j,k,it)*E1(p,i,j,k);
	  }}
	  for (j=0; j<Ntan(p,0)-1; j++){
	    for (k=0; k<Ntan(p,1); k++){ 		
	      E2(p,i,j,k)     = Ae(p,i,0)*E2(p,i,j,k) + Ae(p,i,1)*(Hx(j,k,it+1)-Hx(j,k,it));
	      DiffEy(j,k,it) += Cey(j,k,it)*E2(p,i,j,k);
	  }}
	}
}



void updateDiffHpml(Grid *g){   
	
	int i, j, k, p, it;  // three spatial indices, PML index, translated i index

	p=0; // left x-pml
	for (i=0; i<Npml(p); i++){  
	  for (j=0; j<Ntan(p,0)-1; j++){
	    for (k=0; k<Ntan(p,1); k++){ 	
	      H1(p,i,j,k)    = Ah(p,i,0)*H1(p,i,j,k) + Ah(p,i,1)*(Ez(i+1,j+1,k)-Ez(i,j+1,k));
	      DiffHy(i,j,k) -= Chy(i,j,k)*H1(p,i,j,k); 
	  }} 
	  for (j=0; j<Ntan(p,0); j++){
	    for (k=0; k<Ntan(p,1)-1; k++){ 
	      H2(p,i,j,k)    = Ah(p,i,0)*H2(p,i,j,k) + Ah(p,i,1)*(Ey(i+1,j,k+1)-Ey(i,j,k+1));  
	      DiffHz(i,j,k) += Chz(i,j,k)*H2(p,i,j,k);	
          }}
	}

	p=1;  // right x-pml	
	for (i=0; i<Npml(p); i++){  
	  it = Nx-Npml(p)+i;  
	  for (j=0; j<Ntan(p,0)-1; j++){
	    for (k=0; k<Ntan(p,1); k++){ 
	      H1(p,i,j,k)     = Ah(p,i,0)*H1(p,i,j,k) + Ah(p,i,1)*(Ez(it+1,j+1,k)-Ez(it,j+1,k));
	      DiffHy(it,j,k) -= Chy(it,j,k)*H1(p,i,j,k); 
	  }} 
	  for (j=0; j<Ntan(p,0); j++){
	    for (k=0; k<Ntan(p,1)-1; k++){ 
	      H2(p,i,j,k)     = Ah(p,i,0)*H2(p,i,j,k) + Ah(p,i,1)*(Ey(it+1,j,k+1)-Ey(it,j,k+1));  
	      DiffHz(it,j,k) += Chz(it,j,k)*H2(p,i,j,k);	
          }}
	}

	p=2; // left y-pml
	for (i=0; i<Npml(p); i++){  
	  for (j=0; j<Ntan(p,0)-1; j++){
	    for (k=0; k<Ntan(p,1); k++){ 	
	      H1(p,i,j,k)    = Ah(p,i,0)*H1(p,i,j,k) + Ah(p,i,1)*(Ex(k,i+1,j+1)-Ex(k,i,j+1));
	      DiffHz(k,i,j) -= Chz(k,i,j)*H1(p,i,j,k); 
	  }} 
	  for (j=0; j<Ntan(p,0); j++){
	    for (k=0; k<Ntan(p,1)-1; k++){ 
	      H2(p,i,j,k)    = Ah(p,i,0)*H2(p,i,j,k) + Ah(p,i,1)*(Ez(k+1,i+1,j)-Ez(k+1,i,j));  
	      DiffHx(k,i,j) += Chx(k,i,j)*H2(p,i,j,k);	
          }}
	}

	p=3; // right y-pml
	for (i=0; i<Npml(p); i++){  
	  it = Ny-Npml(p)+i;
	  for (j=0; j<Ntan(p,0)-1; j++){
	    for (k=0; k<Ntan(p,1); k++){ 	 
	      H1(p,i,j,k)     = Ah(p,i,0)*H1(p,i,j,k) + Ah(p,i,1)*(Ex(k,it+1,j+1)-Ex(k,it,j+1));
	      DiffHz(k,it,j) -= Chz(k,it,j)*H1(p,i,j,k); 
	  }} 
	  for (j=0; j<Ntan(p,0); j++){
	    for (k=0; k<Ntan(p,1)-1; k++){ 
	      H2(p,i,j,k)     = Ah(p,i,0)*H2(p,i,j,k) + Ah(p,i,1)*(Ez(k+1,it+1,j)-Ez(k+1,it,j));  
	      DiffHx(k,it,j) += Chx(k,it,j)*H2(p,i,j,k);	
          }}
	}

	p=4; // left z-pml
	for (i=0; i<Npml(p); i++){  
	  for (j=0; j<Ntan(p,0)-1; j++){
	    for (k=0; k<Ntan(p,1); k++){ 	
	      H1(p,i,j,k)    = Ah(p,i,0)*H1(p,i,j,k) + Ah(p,i,1)*(Ey(j+1,k,i+1)-Ey(j+1,k,i));
	      DiffHx(j,k,i) -= Chx(j,k,i)*H1(p,i,j,k); 
	  }} 
	  for (j=0; j<Ntan(p,0); j++){
	    for (k=0; k<Ntan(p,1)-1; k++){ 
	      H2(p,i,j,k)    = Ah(p,i,0)*H2(p,i,j,k) + Ah(p,i,1)*(Ex(j,k+1,i+1)-Ex(j,k+1,i));  
	      DiffHy(j,k,i) += Chy(j,k,i)*H2(p,i,j,k);	
          }}
	}

	p=5; // right z-pml
	for (i=0; i<Npml(p); i++){  
	  it = Nz-Npml(p)+i;
	  for (j=0; j<Ntan(p,0)-1; j++){
	    for (k=0; k<Ntan(p,1); k++){ 	
	      H1(p,i,j,k)     = Ah(p,i,0)*H1(p,i,j,k) + Ah(p,i,1)*(Ey(j+1,k,it+1)-Ey(j+1,k,it));
	      DiffHx(j,k,it) -= Chx(j,k,it)*H1(p,i,j,k); 
	  }} 
	  for (j=0; j<Ntan(p,0); j++){
	    for (k=0; k<Ntan(p,1)-1; k++){ 
	      H2(p,i,j,k)     = Ah(p,i,0)*H2(p,i,j,k) + Ah(p,i,1)*(Ex(j,k+1,it+1)-Ex(j,k+1,it));  
	      DiffHy(j,k,it) += Chy(j,k,it)*H2(p,i,j,k);	
          }}
	}
}


void freeMemoryPml(){
	for(int p=5; p>=06; p--){
	  free(pmls[p].ah); free(pmls[p].ae);  
	  free(pmls[p].h2); free(pmls[p].h1); 
	  free(pmls[p].e2); free(pmls[p].e1);
	}
}



