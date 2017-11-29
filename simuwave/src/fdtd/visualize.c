
#include "fdtd3d.h"
#include <math.h>  


// static variables preserve their value in between invocations and are only seen by this file
static int type, start, stride;     
static double q; 
static PyObject *vfun, *plist;


void initVisualization(Grid *g, PyObject *ipt){
	
	PyObject *py_vis = PyObject_GetAttrString( ipt, "vis" );

	// field type: {0,1,2,3,4} = {nothing,Ex,Ey,Ez,E} 
	PyObject *py_type = PyObject_GetAttrString( py_vis, "type" );	
	type = PyInt_AS_LONG(py_type);	
	Py_DECREF(py_type);

	if (type!=0){	

		// time index of the first scene
		PyObject *py_start = PyObject_GetAttrString( py_vis, "start" );
		start = PyInt_AS_LONG(py_start);
		Py_DECREF(py_start);
		
		// number of iterations between two scenes
		PyObject *py_stride = PyObject_GetAttrString( py_vis, "stride" );
		stride = PyInt_AS_LONG(py_stride);
		Py_DECREF(py_stride);

		// rescaling factor (python-float = C-double)
		PyObject *py_scale = PyObject_GetAttrString( py_vis, "scale" );
                q =  Z0/4.0/PyFloat_AS_DOUBLE(py_scale);
		Py_DECREF(py_scale);                                   	 
		
	        // initialize Python list used in the visualize function below    
		plist = PyList_New( (Nx-Npml(0)-Npml(1))* \
                                    (Ny-Npml(2)-Npml(3))* \
                                    (Nz-Npml(4)-Npml(5)) );
	
		// append the directory where 'visualize.py' is located to sys.path
		char *argv[] = {"visualize.py"};        // array of strings 
		int argc = 1;                           // number of strings in the array on the previous line 		
		PySys_SetArgvEx(argc,argv,1);                               
		
		// import and apply the python module and its functions
		PyObject *mod, *ifun, *mg, *arg;		
		mod   = PyImport_ImportModule("visualize");                 	
		ifun  = PyObject_GetAttrString(mod, "initVisualization");    
		vfun  = PyObject_GetAttrString(mod, "visualize");
		mg    = PyObject_GetAttrString(ipt, "maingrid");							
		arg   = Py_BuildValue("(O)", mg);       // construct argument (a Python tuple) --> new reference		
		PyEval_CallObject(ifun, arg);                                
		Py_DECREF(arg);                                           			 
		Py_DECREF(mg);
		Py_DECREF(ifun);
		Py_DECREF(mod); 

	}

	Py_DECREF(py_vis);

}



void visualize(Grid *g){
	
	if (type!=0 && It>=start && (It-start)%stride==0){
		
		PyObject *arg;
		Py_ssize_t index = 0;  
			
		int i, j, k;
                                  		
                switch (type){
	          case 1:
	            for(i=Npml(0); i<Nx-Npml(1); i++){ 
		      for(j=Npml(2); j<Ny-Npml(3); j++){
		        for(k=Npml(4); k<Nz-Npml(5); k++){				  
			  my_SetItem( plist, index, q*(Ex(i,j,k)+Ex(i,j+1,k)+Ex(i,j,k+1)+Ex(i,j+1,k+1))/Dx1(i) );    
			  index++;  
                    }}}  
		    break; 
		  case 2:
	            for(i=Npml(0); i<Nx-Npml(1); i++){ 
		      for(j=Npml(2); j<Ny-Npml(3); j++){
		        for(k=Npml(4); k<Nz-Npml(5); k++){				  
			  my_SetItem( plist, index, q*(Ey(i,j,k)+Ey(i+1,j,k)+Ey(i,j,k+1)+Ey(i+1,j,k+1))/Dy1(j) );    
			  index++;  
                    }}}  
		    break; 
                  case 3:
	            for(i=Npml(0); i<Nx-Npml(1); i++){ 
		      for(j=Npml(2); j<Ny-Npml(3); j++){
		        for(k=Npml(4); k<Nz-Npml(5); k++){			  
			  my_SetItem( plist, index, q*(Ez(i,j,k)+Ez(i+1,j,k)+Ez(i,j+1,k)+Ez(i+1,j+1,k))/Dz1(k) );   
			  index++;  
                    }}}  
		    break; 
                  case 4:
	            for(i=Npml(0); i<Nx-Npml(1); i++){ 
		      for(j=Npml(2); j<Ny-Npml(3); j++){
		        for(k=Npml(4); k<Nz-Npml(5); k++){						  
			  my_SetItem( plist, index, q*sqrt( pow( (Ex(i,j,k)+Ex(i,j+1,k)+Ex(i,j,k+1)+Ex(i,j+1,k+1))/Dx1(i) , 2) + \
					                    pow( (Ey(i,j,k)+Ey(i+1,j,k)+Ey(i,j,k+1)+Ey(i+1,j,k+1))/Dy1(j) , 2) + \
					                    pow( (Ez(i,j,k)+Ez(i+1,j,k)+Ez(i,j+1,k)+Ez(i+1,j+1,k))/Dz1(k) , 2) ) );
			  index++;  
                    }}}  
		}
	
		arg = Py_BuildValue("(O)", plist);         	
		PyEval_CallObject(vfun, arg);
		Py_DECREF(arg);  		
 		
	}
}




void freeMemoryVisualization(){

	if (type!=0){ 
	  Py_DECREF(plist); 
	  Py_DECREF(vfun); 
        }
	
}



/* 
COMMENTS:
(1) Python is embedded inside the C-funtion using the Python-C API, such that the Mayavi 3D renderer can 
      be used inside the C timestepping loop (avoiding the need for an expensive interpreted Python for-loop).
(2) PyList_SET_ITEM is faster than PyList_SetItem, but more prone to errors (memory leaks). 
      Both have stolen references and do not require the use of Py_DECREF.
(3) Actually, the electric field components should be interpolated to the centers of the dual 
      Yee cells in order to preserve second-order accuracy. 
(4) The field data should be mapped to a grid with uniform step equal to the largest step 
      occurring in the nonuniform grid.
(5) In Py_BuildValue: "O" increments the reference count, whereas "N" does not. 
*/




/*
TERMINAL:
sudo find / -name "Python.h"   -->    /usr/include/python2.7/Python.h
gcc <<c-files>> -I/usr/include/python2.7 -lpython2.7

Explanation:
The argument -I<<path>> adds a new path to the gcc search path. 
The argument -lpython2.7 is used to link the Python library with the executable.
*/

