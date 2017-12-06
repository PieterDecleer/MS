#include "fdtd3d.h"
#include <math.h>

void initSchrodgrid(Grid *g, PyObject *ipt){
	printf("Initialization of the Schrodinger grid started...\n");
	//definition of the reduced mass
	static double mr = (0.023*Me);
	
	// three spatial indices and one subgrid index
	int i,j,k,s;

	// number of Schrodinger grids
	Nsr = 1;	// If multiple grids are used with subgridding, this has to be obtained from the GUI
	
	for (s=0; s<Nsr; s++) {

		// number of cells
		Nxs = 50;
		Nys = 50;
		Nzs = 50;
		int xc = Nxs/2;							// center of harmonic oscillator in x direction
		int yc = Nys/2;							// center of harmonic oscillator in y direction
		int zc = Nys/2;							// center of harmonic oscillator in z direction

		// redefine timestep to normal time, not Minkowski time
		Dts = Dtau/c0;
		
		// cell dimensions
		my_alloc( (*g).srg[s].dxs, Nxs, double);
		my_alloc( (*g).srg[s].dys, Nys, double);
		my_alloc( (*g).srg[s].dzs, Nzs, double);
		for (i=0; i<Nxs; i++) Dxs(i) = Dx1(0)/Nxs;
		for (j=0; j<Nys; j++) Dys(j) = Dy1(0)/Nys;
		for (k=0; k<Nzs; k++) Dzs(k) = Dz1(0)/Nzs;

		printf("Cell size in center = %e x %e x %e\n",Dxs(xc),Dys(yc),Dzs(zc));
 		
		Dts = 0.1/((Hbar/mr)*(12/(Dx1(xc)*Dx1(xc))));

		// dynamically allocate memoryfor the field variables and the update coefficients
		my_alloc( (*g).srg[s].pr     , Nxs*Nys*Nzs     , double );	// Ask for indices
		my_alloc( (*g).srg[s].prv    , Nxs*Nys*Nzs     , double );
		my_alloc( (*g).srg[s].prvv   , Nxs*Nys*Nzs     , double );
		my_alloc( (*g).srg[s].pi     , Nxs*Nys*Nzs     , double );
		my_alloc( (*g).srg[s].piv    , Nxs*Nys*Nzs     , double );
		my_alloc( (*g).srg[s].pivv   , Nxs*Nys*Nzs     , double );
		my_alloc( (*g).srg[s].lpx    , Nxs     , double );
		my_alloc( (*g).srg[s].lpy    , Nys     , double );
		my_alloc( (*g).srg[s].lpz    , Nzs     , double );
		my_alloc( (*g).srg[s].ppx    , Nxs     , double );
		my_alloc( (*g).srg[s].ppy    , Nys     , double );
		my_alloc( (*g).srg[s].ppz    , Nzs     , double );
		my_alloc( (*g).srg[s].vs     , Nxs*Nys*Nzs     , double );
		
		// set laplacian and potential update coefficients
		for (i=1; i<Nxs; i++){
			Lpx(i) = Hbar*Dts/(mr*Dxs(xc)*Dxs(xc));		// for a non uniform grid check what has to be used!!!
			Ppx(i) = 2*Dts/Hbar;
		}
		for (j=1; j<Nys; j++){
			Lpy(j) = Hbar*Dts/(mr*Dys(yc)*Dys(yc));		// for a non uniform grid check what has to be used!!!
			Ppy(i) = 2*Dts/Hbar;
		}
		for (k=1; k<Nzs; k++){
			Lpz(k) = Hbar*Dts/(mr*Dzs(zc)*Dzs(zc));		// for a non uniform grid check what has to be used!!!
			Ppz(k) = 2*Dts/Hbar;
		}

		// set static potential strength and initial wavefunction
		static double omega = (2*M_PI*(c0))/(950e-9);		// static potential angular frequency

		
		double xdist(int pos, int xc){
			if (pos>xc)
				return xdist(pos-1, xc)+Dxs(pos-1);
			else if (pos<xc)
				return xdist(pos+1, xc)+Dxs(pos);
			else return 0.0;
		}
		double ydist(int pos, int yc){
			if (pos>yc)
				return ydist(pos-1, yc)+Dys(pos-1);
			else if (pos<yc)
				return ydist(pos+1, yc)+Dys(pos);
			else return 0.0;
		}
		double zdist(int pos, int zc){
			if (pos>zc)
				return zdist(pos-1, zc)+Dzs(pos-1);
			else if (pos<zc)
				return zdist(pos+1, zc)+Dzs(pos);
			else return 0.0;
		}
		double coeff = pow(mr*omega/(M_PI*Hbar),0.75);
		double coeff2 = mr*omega/(2.0*Hbar);
		for (i=0; i<Nxs; i++){
			for (j=0; j<Nys; j++){
				for (k=0; k<Nzs; k++){
					double distance_squared = pow((i-xc)*Dxs(xc), 2.0) + pow((j-yc)*Dys(yc), 2.0) + pow((k-zc)*Dzs(zc), 2.0);
					Vs(i,j,k) = 0.5*mr*pow(omega, 2.0)*distance_squared;
					Prvv(i,j,k) = coeff*exp(-coeff2*distance_squared);
					Pivv(i,j,k) = coeff*exp(-coeff2*distance_squared)*0.0;
					Prv(i,j,k) = coeff*exp(-coeff2*distance_squared)*cos(-1.5*omega*Dts);
					Piv(i,j,k) = coeff*exp(-coeff2*distance_squared)*sin(-1.5*omega*Dts);
					
					
					 		
		}}}

		
		FILE *f = fopen("file.txt","w");
		if (f == NULL){
			printf("Error opening file!\n");
			exit(1);
		}			
		for (i=0; i<Nxs; i++){
			for (j =0; j<Nys; j++){
				fprintf(f,"%e\t", Prvv(i,j,zc));
		} fprintf(f,"\n");}		
		fclose(f);
		printf("coeff = %e\ncoeff2 = %e\n", coeff, coeff2);
	}
	return;
}

void updateElectron(Grid *g){
	
    int i,j,k,s;
	for (s=0;s<Nsr;s++){
		for (i = 1; i < Nxs - 1; i++){
		    for (j = 1; j < Nys - 1; j++){
		        for (k = 1; k < Nzs - 1; k++){
		            Pr(i,j,k) = Prvv(i,j,k)     // update Pr to time n+1/2
		                        - Lpx(i) * (Piv(i + 1,j,k) - 2 * Piv(i,j,k) + Piv(i - 1,j,k))
		                        - Lpy(j) * (Piv(i,j + 1,k) - 2 * Piv(i,j,k) + Piv(i,j - 1,k))
		                        - Lpz(k) * (Piv(i,j,k + 1) - 2 * Piv(i,j,k) + Piv(i,j,k - 1))
		                        + Ppx(i) * Vs(i,j,k) * Piv(i,j,k);
		        }
		    }
		}
		for (i = 1; i < Nxs - 1; i++){
		    for (j = 1; j < Nys - 1; j++){
		        for (k = 1; k < Nzs - 1; k++){
		            Pi(i,j,k) = Pivv(i,j,k)     // update Pi to time n+1/2
		                        + Lpx(i) * (Prv(i + 1,j,k) - 2 * Prv(i,j,k) + Prv(i - 1,j,k))
		                        + Lpy(j) * (Prv(i,j + 1,k) - 2 * Prv(i,j,k) + Prv(i,j - 1,k))
		                        + Lpz(k) * (Prv(i,j,k + 1) - 2 * Prv(i,j,k) + Prv(i,j,k - 1))
		                        - Ppx(i) * Vs(i,j,k) * Prv(i,j,k);
		        }
		    }
		}
		for (i = 1; i < Nxs - 1; i++){
		    for (j = 1; j < Nys - 1; j++){
		        for (k = 1; k < Nzs - 1; k++){
		            Prvv(i,j,k) = Prv(i,j,k);   // update Prvv to time n-1/2
		        }
		    }
		}
		for (i = 1; i < Nxs - 1; i++){
		    for (j = 1; j < Nys - 1; j++){
		        for (k = 1; k < Nzs - 1; k++){
		            Pivv(i,j,k) = Piv(i,j,k);   // update Pivv to time n-1/2
		        }
		    }
		}

		for (i = 1; i < Nxs - 1; i++){
		    for (j = 1; j < Nys - 1; j++){
		        for (k = 1; k < Nzs - 1; k++){
		            Prv(i,j,k) = Pr(i,j,k);     // update Prv to time n+1/2
		        }
		    }
		}
		for (i = 1; i < Nxs - 1; i++){
		    for (j = 1; j < Nys - 1; j++){
		        for (k = 1; k < Nzs - 1; k++){
		            Piv(i,j,k) = Pi(i,j,k);     // update Piv to time n+1/2
		        }
		    }
		}


	}
	
    return;
}


int frame = 0;

void electronSnapshot(Grid *g){
	int i,j,s;	
	FILE *snapshot;
	char f[100];

	if (It % 100 == 0){
		printf("snapshot taken\n");
		for (s=0;s<Nsr;s++){
			sprintf(f, "snapshot_%d.txt", frame++);
			snapshot = fopen(f, "w");
			int zc = Nys/2;			
			for (i=0; i<Nxs; i++){
				for (j =0; j<Nys; j++){
					fprintf(snapshot,"%e\t", pow(Pr(i,j,zc),2.0)+pow(Pi(i,j,zc),2.0));
			} fprintf(snapshot,"\n");}		
			fclose(snapshot);
		}
	}
}

void freeMemorySchrodgrid(Grid *g){
	int s;
	for (s=Nsr-1;s>=0;s--){ 
		free((*g).srg[s].vs);
		free((*g).srg[s].ppz); free((*g).srg[s].ppy); free((*g).srg[s].ppx);
		free((*g).srg[s].lpz); free((*g).srg[s].lpy); free((*g).srg[s].lpx);
		free((*g).srg[s].pivv); free((*g).srg[s].piv); free((*g).srg[s].pi);
		free((*g).srg[s].prvv); free((*g).srg[s].prv); free((*g).srg[s].pr);
		free((*g).srg[s].dzs); free((*g).srg[s].dys); free((*g).srg[s].dxs);
	}
}




