
#include "fdtd3d.h"                                                                


PyObject* run_fdtd(PyObject *ipt){

	printf("Initializing ...\n"); 

	Grid *g; my_alloc(g,1,Grid);
	
	initGrid(g,ipt);
	initAdhie1(g,ipt);
	initPml(g,ipt);
	initAdhie2(g); 
	initSource(g,ipt);	 	
	initSensor(g,ipt);	
	initVisualization(g,ipt);  	
	
	printf("Timestepping ...\n"); 

	for (It=0; It<Nt; It++){                         
	
		printf("%i/%i\n",It,Nt); 		
		
		updateDiffH(g);
		updateDiffHpml(g);		
		updateDiffHadhie(g); 		
		updateH(g); 
		updateHhard(g); 					

		//updateSubgrids(g);

		updateDiffE(g);
		//fine2coarse(g);
		updateDiffEcurrent(g);	
		updateDiffEpml(g);
		updateDiffEadhie(g);		
		updateE(g);	
		updateEhard(g);
 	
		//coarse2fine(g);	
		
		updateSensor(g); 
		visualize(g);	
		
	}

	printf("Free memory ...\n"); 

	freeMemoryVisualization();
	PyObject *out = sensorC2Py(g);
	freeMemorySource(); 
	freeMemoryPml();
	freeMemoryAdhie();     
	freeMemoryGrid(g);   

	return out; 
}


/* 

Execution in terminal for standalone C-code:
gcc -Wall -O -c <<list of all c-files>>
gcc <<list of all o-files>> -lm -o fdtd3d.exe 
./fdtd3d.exe

Explanation: 
  Wall: enable all warnings
  O   : optimization for code size and execution time (O2 and O3 optimize even more at the cost of a higher compilation time) 
  c   : compile into object file (o-file) without linking
  lm  : link to the math library
  o   : write the build output to a specified output file (here fdtd3d.exe)	

Including Python visualization: 
gcc -Wall -O -c fdtd3d.c grid.c pml.c source.c visualize.c -I/usr/include/python2.7
gcc fdtd3d.o grid.o pml.o source.o visualize.o -lm -lpython2.7 -o fdtd3d.exe
./fdtd3d.exe

*/
