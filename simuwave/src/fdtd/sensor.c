
#include "fdtd3d.h"
#include <math.h>

static double *sensor;                           
static int Ns, *type, *i_1, *j_1, *k_1, *i_2, *j_2, *k_2;


void initSensor(Grid *g, PyObject *ipt){
            
	PyObject *py_mg     = PyObject_GetAttrString( ipt, "maingrid" );
	PyObject *py_sensor = PyObject_GetAttrString( py_mg, "sensor" );            
	
	// number of sensors
	Ns = (int) PyList_Size(py_sensor);   
	
	// dynamic memory allocation
	my_alloc(sensor, Ns*Nt, double); 
	my_alloc(type, Ns, int);	
	my_alloc(i_1 , Ns, int);
	my_alloc(i_2 , Ns, int);
	my_alloc(j_1 , Ns, int);
	my_alloc(j_2 , Ns, int);
	my_alloc(k_1 , Ns, int);
	my_alloc(k_2 , Ns, int);
 
	for(int s=0; s<Ns; s++){
	  
    	  // sensor type
          PyObject *py_type = PyObject_GetAttrString(PyList_GET_ITEM (py_sensor, (Py_ssize_t) s), "type" );
	  type[s] = PyInt_AS_LONG(py_type);
 	  Py_DECREF(py_type);
	  
          // sensor location
          PyObject *py_ind  = PyObject_GetAttrString(PyList_GET_ITEM( py_sensor, (Py_ssize_t) s), "ind" );
	  i_1[s]  = PyInt_AS_LONG( PyList_GET_ITEM( py_ind, 0) );
	  j_1[s]  = PyInt_AS_LONG( PyList_GET_ITEM( py_ind, 1) );
	  k_1[s]  = PyInt_AS_LONG( PyList_GET_ITEM( py_ind, 2) );
	  i_2[s]  = PyInt_AS_LONG( PyList_GET_ITEM( py_ind, 3) );
	  j_2[s]  = PyInt_AS_LONG( PyList_GET_ITEM( py_ind, 4) );
	  k_2[s]  = PyInt_AS_LONG( PyList_GET_ITEM( py_ind, 5) );
	  Py_DECREF(py_ind);

	}

	Py_DECREF(py_sensor);
	Py_DECREF(py_mg);  

}



void updateSensor(Grid *g){

	int i,j,k;

	for(int s=0; s<Ns; s++){
		  
  	  switch (type[s]){
	     
          // magnetic field [A/m]	
            case 0: sensor[Ns*It+s] = Hx(i_1[s],j_1[s],k_1[s])/Dx2(i_1[s]); break;         
            case 1: sensor[Ns*It+s] = Hy(i_1[s],j_1[s],k_1[s])/Dy2(j_1[s]); break;
            case 2: sensor[Ns*It+s] = Hz(i_1[s],j_1[s],k_1[s])/Dz2(k_1[s]); break;

          // electric field [V/m]
	    case 3: sensor[Ns*It+s] = Z0*Ex(i_1[s],j_1[s]+1,k_1[s]+1)/Dx1(i_1[s]); break;  
	    case 4: sensor[Ns*It+s] = Z0*Ey(i_1[s]+1,j_1[s],k_1[s]+1)/Dy1(j_1[s]); break;
            case 5: sensor[Ns*It+s] = Z0*Ez(i_1[s]+1,j_1[s]+1,k_1[s])/Dz1(k_1[s]); break;

          // quasi-static voltage [V]
	    case 6:                                                            
	      for(i=i_1[s]; i<=i_2[s]; i++){
                  sensor[Ns*It+s] += Z0*Ex(i,j_1[s]+1,k_1[s]+1); 
              } break;
            case 7:
	      for(j=j_1[s]; j<=j_2[s]; j++){
                  sensor[Ns*It+s] += Z0*Ey(i_1[s]+1,j,k_1[s]+1); 
              } break;
	    case 8:
	      for(k=k_1[s]; k<=k_2[s]; k++){
                  sensor[Ns*It+s] += Z0*Ez(i_1[s]+1,j_1[s]+1,k); 
              } break; 
 
          // current [A] (factor Z0??)
	    case 9:                                                             
	      for(j=j_1[s]; j<=j_2[s]; j++){
	        for(k=k_1[s]; k<=k_2[s]; k++){
                  sensor[Ns*It+s] -= Sex(i_1[s],j,k)/Cex(i_1[s],j,k)*Ex(i_1[s],j+1,k+1); 
              }} break;
            case 10:
	      for(i=i_1[s]; i<=i_2[s]; i++){
	        for(k=k_1[s]; k<=k_2[s]; k++){
                  sensor[Ns*It+s] -= Sey(i,j_1[s],k)/Cey(i,j_1[s],k)*Ey(i+1,j_1[s],k+1); 
              }} break;
	    case 11:	
	      for(i=i_1[s]; i<=i_2[s]; i++){
	        for(j=j_1[s]; j<=j_2[s]; j++){
                  sensor[Ns*It+s] -= Sez(i,j,k_1[s])/Cez(i,j,k_1[s])*Ez(i+1,j+1,k_1[s]); 
              }} break; 

           // quasi-static current [A] 
	     case 12:                                                          
	       for(j=j_1[s]; j<=j_2[s]+1; j++){ sensor[Ns*It+s] += (Hy(i_1[s],j,k_1[s])-Hy(i_1[s],j,k_2[s]+1)); }
	       for(k=k_1[s]; k<=k_2[s]+1; k++){ sensor[Ns*It+s] += (Hz(i_1[s],j_2[s]+1,k)-Hz(i_1[s],j_1[s],k)); } break; 	
	     case 13:	
	       for(i=i_1[s]; i<=i_2[s]+1; i++){ sensor[Ns*It+s] += (Hx(i,j_1[s],k_2[s]+1)-Hx(i,j_1[s],k_1[s])); }
	       for(k=k_1[s]; k<=k_2[s]+1; k++){ sensor[Ns*It+s] += (Hz(i_1[s],j_1[s],k)-Hz(i_2[s]+1,j_1[s],k)); } break;
	     case 14:
	       for(i=i_1[s]; i<=i_2[s]+1; i++){ sensor[Ns*It+s] += (Hx(i,j_1[s],k_1[s])-Hx(i,j_2[s]+1,k_1[s])); }	
	       for(j=j_1[s]; j<=j_2[s]+1; j++){ sensor[Ns*It+s] += (Hy(i_2[s]+1,j,k_1[s])-Hy(i_1[s],j,k_1[s])); } 
	  }
        }
}



PyObject* sensorC2Py(Grid *g){	
	Py_ssize_t N = Ns*Nt; 	
	PyObject *plist = PyList_New(N); 
	for(Py_ssize_t i=0; i<N; i++){ my_SetItem(plist, i, sensor[i]); }
	free(sensor); 
        return plist; 
}




/* COMMENTS:
(1) EMPro: Current densities are determined by multiplying the conductivity of the material at the specified location 
      by the electric field in the given direction. When a PEC material is present, the current density will be computed 
      by the loop of magnetic fields surrounding that cell edge. 
(2) Current integral coefficients should be precomputed and dynamically stored.  
(3) For electric fields located on the material boundaries, the contribution to the current surface integral should be
      divided by two, which is automatically taken into account by the averaged sigma. 
*/


