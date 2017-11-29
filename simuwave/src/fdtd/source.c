
#include "fdtd3d.h"


static double *source;                           
static int type, i_1, j_1, k_1, i_2, j_2, k_2;


void initSource(Grid *g, PyObject *ipt){
        
        // dynamic memory allocation
	PyObject *py_wf = PyObject_GetAttrString( ipt, "wf" );
	PyObject *py_ts = PyObject_GetAttrString( py_wf, "timesignal" );        
	my_alloc(source, Nt, double); 
	for (int i=0; i<Nt; i++){ source[i] = PyFloat_AS_DOUBLE( PyList_GET_ITEM( py_ts, (Py_ssize_t) i) ); } // definition of the source
	Py_DECREF(py_ts);
	Py_DECREF(py_wf);

	PyObject *py_mg     = PyObject_GetAttrString( ipt, "maingrid" );
	PyObject *py_source = PyObject_GetAttrString( py_mg, "source" );	
	
	if (py_source != Py_None){	
	 
	  // source type	
          PyObject *py_type = PyObject_GetAttrString( py_source, "type");
	  type = PyInt_AS_LONG(py_type);	
	  Py_DECREF(py_type);
           
          // source position
	  PyObject *py_ind = PyObject_GetAttrString(py_source, "ind");             
	  i_1 = PyInt_AS_LONG( PyList_GET_ITEM( py_ind, 0) );
	  j_1 = PyInt_AS_LONG( PyList_GET_ITEM( py_ind, 1) );
	  k_1 = PyInt_AS_LONG( PyList_GET_ITEM( py_ind, 2) );
	  i_2 = PyInt_AS_LONG( PyList_GET_ITEM( py_ind, 3) );
	  j_2 = PyInt_AS_LONG( PyList_GET_ITEM( py_ind, 4) );
	  k_2 = PyInt_AS_LONG( PyList_GET_ITEM( py_ind, 5) );
	  Py_DECREF(py_ind); 	
	  
	}
 
	Py_DECREF(py_source);
	Py_DECREF(py_mg);

	if (type>=3 && type<=5) { for (int i=0; i<Nt; i++){ source[i] /= Z0; } }

}



void updateHhard(Grid *g){

	int i,j,k;

	switch (type){
          case 0:	
            for(i=i_1; i<=i_2; i++){
	      for(j=j_1; j<=j_2; j++){
	        for(k=k_1; k<=k_2; k++){		
		  Hx(i,j,k) = Dx2(i)*source[It];
            }}} 
            break; 
          case 1: 
            for(i=i_1; i<=i_2; i++){
	      for(j=j_1; j<=j_2; j++){
	        for(k=k_1; k<=k_2; k++){		
		  Hy(i,j,k) = Dy2(j)*source[It];
            }}} 
            break; 
          case 2: 
            for(i=i_1; i<=i_2; i++){
	      for(j=j_1; j<=j_2; j++){
	        for(k=k_1; k<=k_2; k++){		
		  Hz(i,j,k) = Dz2(k)*source[It];
            }}} 	  		            
	} 
}



void updateEhard(Grid *g){

	int i,j,k;

	switch (type){
	  case 3: 
	    for(i=i_1; i<=i_2; i++){
	      for(j=j_1; j<=j_2; j++){
	        for(k=k_1; k<=k_2; k++){		
		  Ex(i,j+1,k+1) = Dx1(i)*source[It]; 
            }}} 
            break;	 
	  case 4: 
	    for(i=i_1; i<=i_2; i++){
	      for(j=j_1; j<=j_2; j++){
	        for(k=k_1; k<=k_2; k++){		
		  Ey(i+1,j,k+1) = Dy1(j)*source[It]; 
            }}} 
            break; 	
          case 5: 
            for(i=i_1; i<=i_2; i++){
	      for(j=j_1; j<=j_2; j++){
	        for(k=k_1; k<=k_2; k++){		
		  Ez(i+1,j+1,k) = Dz1(k)*source[It];
            }}} 
	 }
}



void updateDiffEcurrent(Grid *g){       // current [A], not current density [A/m²]       
		
	int i,j,k; 

	switch (type){
	  case 6:  
	    for(i=i_1; i<=i_2; i++){
	      for(j=j_1; j<=j_2; j++){
	        for(k=k_1; k<=k_2; k++){		
		  DiffEx(i,j,k) -= Cex(i,j,k)*source[It];
            }}} 
	    break; 
	  case 7:  
            for(i=i_1; i<=i_2; i++){
	      for(j=j_1; j<=j_2; j++){
	        for(k=k_1; k<=k_2; k++){		
		  DiffEy(i,j,k) -= Cey(i,j,k)*source[It];
            }}} 
            break;   	
          case 8: 
	    for(i=i_1; i<=i_2; i++){
	      for(j=j_1; j<=j_2; j++){
	        for(k=k_1; k<=k_2; k++){		
		  DiffEz(i,j,k) -= Cez(i,j,k)*source[It];
            }}}  
	}
}



void freeMemorySource(){
	free(source); 
}


