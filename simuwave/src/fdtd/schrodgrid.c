#include "fdtd3d.h"
#include <math.h>

void initSchrodgrid(Grid *g, PyObject *ipt){
	//definition of the reduced mass
	static double mr = 0.023*Me;
	
	// three spatial indices and one subgrid index
	int i,j,k,s;

	// number of subgrids
	Nsr = 1;	// If multiple grids are used with subgridding, this has to be obtained from the GUI
	
	for (s=0; s<Nsr; s++) {

		// number of cells
		Nxs = Nx;
		Nys = Ny;
		Nzs = Nz;

		// redefine timestep to normal time, not Minkowski time
		Dts = Dtau/c0;

		// cell dimensions
		my_alloc( (*g).srg[s].dxs, Nxs, double);
		my_alloc( (*g).srg[s].dys, Nys, double);
		my_alloc( (*g).srg[s].dzs, Nzs, double);
		for (i=0; i<Nxs; i++) Dxs(i) = Dx1(i);
		for (j=0; j<Nys; j++) Dys(j) = Dy1(j);
		for (k=0; k<Nzs; k++) Dzs(k) = Dz1(k);
 
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
			Lpx(i) = Hbar*Dts/(mr*Dxs(i)*Dxs(i-1));		// for a non uniform grid check what has to be used!!!
			Ppx(i) = 2*Dts/Hbar;
		}
		for (j=1; j<Nys; j++){
			Lpy(j) = Hbar*Dts/(mr*Dys(j)*Dys(j-1));		// for a non uniform grid check what has to be used!!!
			Ppy(i) = 2*Dts/Hbar;
		}
		for (k=1; k<Nzs; k++){
			Lpz(k) = Hbar*Dts/(mr*Dzs(k)*Dzs(k-1));		// for a non uniform grid check what has to be used!!!
			Ppz(k) = 2*Dts/Hbar;
		}

		// set static potential strength and initial wavefunction
		static double omega = 2*3.1415*c0/(950e-9);		// static potential angular frequency
		int xc = Nxs/2;							// center of harmonic oscillator in x direction
		int yc = Nys/2;							// center of harmonic oscillator in y direction
		int zc = Nys/2;							// center of harmonic oscillator in z direction
		
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
		
		for (i=0; i<Nxs; i++){
			for (j=0; j<Nys; j++){
				for (k=0; k<Nzs; k++){
					double distance_squared = pow(xdist(i, xc), 2.0) +
											  pow(ydist(j, yc), 2.0) +
											  pow(zdist(k, zc), 2.0);
					Vs(i,j,k) = 0.5*mr*pow(omega, 2.0)*distance_squared;
					Pr(i,j,k) = pow(mr*omega/(M_PI*Hbar),0.75)*exp(mr*omega*distance_squared/(2.0*Hbar));
					Prvv(i,j,k) = pow(mr*omega/(M_PI*Hbar),0.75)*exp(mr*omega*distance_squared/(2.0*Hbar)); 		
		}}}

	
	char basename[80] = "elProb_init";
	char filename[100];
	FILE *snapshot;
	if (1){
		sprintf(filename, "%s_%d", basename, s);
		snapshot = fopen(filename, "w");	
		for (i=0; i < Nxs; i++){
			for (j=0; i < Nys; j++){
				fprintf(snapshot, "%e\t", pow(Pr(i,j,zc), 2.0));
		} fprintf (snapshot, "\n"); }
		fclose(snapshot);	
	}
	}
	return;
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




