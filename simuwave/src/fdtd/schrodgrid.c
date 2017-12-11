#include "fdtd3d.h"
#include <math.h>
int frame;
void initSchrodgrid(Grid *g, PyObject *ipt){
	printf("Initialization of the Schrodinger grid started...\n");
	// set number of frames back to zero at the initialization of the schrodgrid
	frame = 0;
	//definition of the reduced mass
	static double mr = (0.023*Me);
	
	// three spatial indices and one subgrid index
	int i,j,k,s;

	// number of Schrodinger grids
	Nsr = 1;	// If multiple grids are used with subgridding, this has to be obtained from the GUI
	
	for (s=0; s<Nsr; s++) {

		// number of cells
		int xw = Nx/2;
		int yw = Ny/2;
		int zw = Nz/2;
		Nxs = 50;
		Nys = 50;
		Nzs = 50;
		int xc = Nxs/2;							// center of harmonic oscillator in x direction
		int yc = Nys/2;							// center of harmonic oscillator in y direction
		int zc = Nzs/2;							// center of harmonic oscillator in z direction
		printf("indices of the center point: %d\t%d\t%d\n", xc, yc, zc);
		
		
		// cell dimensions
		my_alloc( (*g).srg[s].dxs, Nxs, double);
		my_alloc( (*g).srg[s].dys, Nys, double);
		my_alloc( (*g).srg[s].dzs, Nzs, double);
		for (i=0; i<Nxs; i++) Dxs(i) = Dx1(0)/Nxs;
		for (j=0; j<Nys; j++) Dys(j) = Dy1(0)/Nys;
		for (k=0; k<Nzs; k++) Dzs(k) = Dz1(0)/Nzs;

		printf("Cell size in center = %e x %e x %e\n",Dxs(xc),Dys(yc),Dzs(zc));
 		
		

		// dynamically allocate memoryfor the field variables and the update coefficients
		my_alloc((*g).srg[s].pr     , Nxs*Nys*Nzs     , double );	// Ask for indices
		my_alloc((*g).srg[s].prv    , Nxs*Nys*Nzs     , double );
		my_alloc((*g).srg[s].prvv   , Nxs*Nys*Nzs     , double );
		my_alloc((*g).srg[s].pi     , Nxs*Nys*Nzs     , double );
		my_alloc((*g).srg[s].piv    , Nxs*Nys*Nzs     , double );
		my_alloc((*g).srg[s].pivv   , Nxs*Nys*Nzs     , double );
		my_alloc((*g).srg[s].lpx    , Nxs     , double );
		my_alloc((*g).srg[s].lpy    , Nys     , double );
		my_alloc((*g).srg[s].lpz    , Nzs     , double );
		my_alloc((*g).srg[s].ppx    , Nxs     , double );
		my_alloc((*g).srg[s].ppy    , Nys     , double );
		my_alloc((*g).srg[s].ppz    , Nzs     , double );
		my_alloc((*g).srg[s].vs     , Nxs*Nys*Nzs	, double );
		my_alloc((*g).srg[s].distance,Nxs*Nys*Nzs	, double );
		

		// set static potential
		static double omega = (2*M_PI*(c0))/(950e-9);		// static potential angular frequency
		double coeff = pow(mr*omega/(M_PI*Hbar),0.75);
		double coeff2 = mr*omega/(2.0*Hbar);
		
		for (i=0; i<Nxs; i++){
			for (j=0; j<Nys; j++){
				for (k=0; k<Nzs; k++){
					Dist(i,j,k) = sqrt(pow((i-xc)*Dxs(i), 2.0) + pow((j-yc)*Dys(j), 2.0) + pow((k-zc)*Dzs(k), 2.0));
					Vs(i,j,k) = 0.5*mr*omega*omega*Dist(i,j,k)*Dist(i,j,k);
					
					
					 		
		}}}
		// Set time step as a function of the potential
		Dts = 0.5/((Hbar/mr)*(2.0/(Dxs(xc)*Dxs(xc))+2.0/(Dys(yc)*Dys(yc))+2.0/(Dys(yc)*Dys(yc)))+Vs(0,0,0)/Hbar);
		printf("Time step for Schrodinger: %es\n", Dts);
		printf("Time step for Maxwell: %es\n", Dtau/c0);
		printf("The minimum time step is: %es\n", min(Dts, Dtau/c0));
		Dtau = min(Dts*c0,Dtau);
		Dts = Dtau/c0;

		// set up ground state
		for (i=0; i<Nxs; i++){
			for (j=0; j<Nys; j++){
				for (k=0; k<Nzs; k++){
					Prvv(i,j,k) = coeff*exp(-coeff2*Dist(i,j,k)*Dist(i,j,k));
					Pivv(i,j,k) = coeff*exp(-coeff2*Dist(i,j,k)*Dist(i,j,k))*0.0;
					Prv(i,j,k) = coeff*exp(-coeff2*Dist(i,j,k)*Dist(i,j,k))*cos(-1.5*omega*Dts);
					Piv(i,j,k) = coeff*exp(-coeff2*Dist(i,j,k)*Dist(i,j,k))*sin(-1.5*omega*Dts);			
					 		
		}}}

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
		
		// write out setup values
		FILE *f = fopen("file.txt","w");
		if (f == NULL){
			printf("Error opening file!\n");
			exit(1);
		}			
		for (i=0; i<Nxs; i++){
			for (j =0; j<Nys; j++){
				fprintf(f,"%e\t", Vs(i,j,zc));
		} fprintf(f,"\n");}		
		fclose(f);
		

		//setup for data out in time
		char filename[100];
		sprintf(filename, "time_evolution_%d.txt", s);
		FILE *time_ev = fopen(filename, "w"); 
		fprintf(time_ev, "%e\t%e\n", Dts, Dxs(xc));	
		fclose(time_ev);
	}
	return;
}

void updateElectron(Grid *g){
	
    int i,j,k,s;
	
		
	for (s=0;s<Nsr;s++){

		int xw, yw, zw;	//position of the quantum well
		xw = Nx/2;
		yw = Ny/2;
		zw = Nz/2;
		if (It % 100 == 0) printf("\t%e\n",Z0*Ex(xw,yw,zw));
		int xc = Nxs/2;							// center of harmonic oscillator in x direction
		int yc = Nys/2;							// center of harmonic oscillator in y direction
		int zc = Nzs/2;
		for (i = 1; i < Nxs - 1; i++){
		    for (j = 1; j < Nys - 1; j++){
		        for (k = 1; k < Nzs - 1; k++){
		            Pr(i,j,k) = Prvv(i,j,k)     // update Pr to time n+1/2
		                        - Lpx(i) * (Piv(i + 1,j,k) - 2 * Piv(i,j,k) + Piv(i - 1,j,k))
		                        - Lpy(j) * (Piv(i,j + 1,k) - 2 * Piv(i,j,k) + Piv(i,j - 1,k))
		                        - Lpz(k) * (Piv(i,j,k + 1) - 2 * Piv(i,j,k) + Piv(i,j,k - 1))
								- (Qe*Dts/Hbar) * (i-xc)*Dxs(i) * Z0 * Ex(xw,yw,zw) * Piv(i,j,k)
								- (Qe*Dts/Hbar) * (j-yc)*Dys(j) * Z0 * Ey(xw,yw,zw) * Piv(i,j,k)
								- (Qe*Dts/Hbar) * (k-zc)*Dzs(k) * Z0 * Ez(xw,yw,zw) * Piv(i,j,k)
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
								+ (Qe*Dts/Hbar) * (i-xc)*Dxs(i) * Z0 * Ex(xw,yw,zw) * Prv(i,j,k)
								+ (Qe*Dts/Hbar) * (j-yc)*Dys(j) * Z0 * Ey(xw,yw,zw) * Prv(i,j,k)
								+ (Qe*Dts/Hbar) * (k-zc)*Dzs(k) * Z0 * Ez(xw,yw,zw) * Prv(i,j,k)
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
void electronTimeevolution(Grid *g){
	int s;	
	FILE *time_ev;
	char f[100];

	if (It % 1 == 0){
		for (s=0;s<Nsr;s++){
			sprintf(f, "time_evolution_%d.txt",s);
			time_ev = fopen(f, "a");
			int xc = Nxs/2;
			int yc = Nys/2;
			int zc = Nzs/2;
			fprintf(time_ev,"%e\t%e\n", Pr(xc,yc,zc),Pi(xc,yc,zc));	
			fclose(time_ev);
		}
	}
}

void freeMemorySchrodgrid(Grid *g){
	int s;
	for (s=Nsr-1;s>=0;s--){ 
		free((*g).srg[s].distance);
		free((*g).srg[s].vs);
		free((*g).srg[s].ppz); free((*g).srg[s].ppy); free((*g).srg[s].ppx);
		free((*g).srg[s].lpz); free((*g).srg[s].lpy); free((*g).srg[s].lpx);
		free((*g).srg[s].pivv); free((*g).srg[s].piv); free((*g).srg[s].pi);
		free((*g).srg[s].prvv); free((*g).srg[s].prv); free((*g).srg[s].pr);
		free((*g).srg[s].dzs); free((*g).srg[s].dys); free((*g).srg[s].dxs);
	}
}




