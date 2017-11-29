#include "fdtd3d.h"

void initSchrodgrid(Grid *g, PyObject *ipt){

	my_alloc( (*g).pr     , Nx*Ny*(Nz)     , double );	// Ask for indices
	my_alloc( (*g).prv    , Nx*Ny*(Nz)     , double );
	my_alloc( (*g).prvv   , Nx*Ny*(Nz)     , double );
	my_alloc( (*g).pi     , Nx*Ny*(Nz)     , double );
	my_alloc( (*g).piv    , Nx*Ny*(Nz)     , double );
	my_alloc( (*g).pivv   , Nx*Ny*(Nz)     , double );
	my_alloc( (*g).lpx    , Nx*Ny*(Nz)     , double );
	my_alloc( (*g).lpy    , Nx*Ny*(Nz)     , double );
	my_alloc( (*g).lpz    , Nx*Ny*(Nz)     , double );
	my_alloc( (*g).ppx    , Nx*Ny*(Nz)     , double );
	my_alloc( (*g).ppy    , Nx*Ny*(Nz)     , double );
	my_alloc( (*g).ppz    , Nx*Ny*(Nz)     , double );
	my_alloc( (*g).vs     , Nx*Ny*Nz       , double );



	return;
}
