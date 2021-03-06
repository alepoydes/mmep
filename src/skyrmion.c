#include "vector.h"
#include "skyrmion.h"
#include "debug.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

// Physical parameters
real* energy_shift_per_atom=NULL; // mean energy of S=(0,0,1) over all positions
real zeeman_shift=NAN;
real anisotropy_shift=NAN;
real* exchange_shift=NULL;
real magnetic_field[3]={0,0,0};
real* nonuniform_field=NULL;
// Structure of crystal lattice
int sizex=0; // Width in unit cells
int sizey=0; // Depth in unit cells
int sizez=0; // Height in unit cells
// Boundary conditions along every axis
int boundary_conditions[3]={-1,-1,-1};
// Tranlation vectors
real translation_vectors[3][3]={{NAN,NAN,NAN},{NAN,NAN,NAN},{NAN,NAN,NAN}};
// Structure of unit cell
int sizeu=0; // Number of atoms in the unit cell
// Position of every atom in the unit cell 
real* atom_positions=NULL;
// Half of number of interactions per unit cell
// Other half is restored by exchange of particles
int sizen=0;
// Relative positions of neighbours for every antom in unit cell
// Every column is of the form
// <x-shift> <y-shift> <z-shift> <s-atom> <d-atom> 
// meaning differences of indices of <d-atom> and <s-atoms> 
// are uqual to (<x-shift>,<y-shift>,<z-shift>)
int* neighbours=NULL;
// Magnetic uniaxial anisotopy K = norm*unit
magnetic_anisotropy_type* magnetic_anisotropy=NULL;
int magnetic_anisotropy_count=0;
// Exchange constant J
real* exchange_constant=NULL;
// dipole interaction constant
real dipole=0;
// Dzyaloshinskii Moriya vector for every pair of atoms
real* dzyaloshinskii_moriya_vector=NULL;
// Spin polarized current
real spin_polarized_current[3]={0,0,0};
real* initial_state=NULL;
int initial_states_count=0;
int* relax_state=NULL;

char* all_active=NULL;
char* active=NULL;
int number_of_active=-1;
int number_of_used=0;
int* positions=NULL;

int dipole_count=0;
int dipole_allocated=0;
int* dipole_idx=NULL; // <source atom> <dest atom> <dest x> <dest y> <dest z>
real* dipole_table=NULL; // <multiplier> <vector x> <vector y> <vector z>

real temperature=0;

int INDEX(int u, int x, int y, int z) {
	return ((u*sizex+x)*sizey+y)*sizez+z;
};

/*
int activate_fast_and_adjacent(const real* grad, int sizep, real threshold) {
	char* data=malloc(SIZE*sizeof(char));
	for(int n=0; n<SIZE; n++) {		
		data[n]=0;
		if(!ISACTIVE(active,n)) continue;
		for(int p=0; p<sizep; p++) {
			const real* d=grad+3*(n+p*SIZE);
			real r=rabs(d[0])+rabs(d[1])+rabs(d[2]);
			if(r>=threshold) { data[n]=1; break; };
		};
	};
	int count=0;
	forall(u,x,y,z) {
		int n=INDEX(u,x,y,z);
		if(!ISACTIVE(all_active, n)) {
			SETPASSIVE(active, n);
			continue;
		};
		char f=data[n];
		if(x>0) f=f||data[INDEX(u,x-1,y,z)];
		if(x<sizex-1) f=f||data[INDEX(u,x+1,y,z)];
		if(y>0) f=f||data[INDEX(u,x,y-1,z)];
		if(y<sizey-1) f=f||data[INDEX(u,x,y+1,z)];
		if(z>0) f=f||data[INDEX(u,x,y,z-1)];
		if(z<sizez-1) f=f||data[INDEX(u,x,y,z+1)];
		if(f) {
			SETACTIVE(active, n);
			count++;
		} else SETPASSIVE(active, n);
	};	
	free(data);
	return count;
};*/

char is_large(const real* grad, int sizep, real threshold) {
	//fprintf(stderr, "is_large(%d, %"RF"g)\n",sizep,RT(threshold));
	for(int p=0; p<sizep; p++) {
		const real* d=grad+p*3*SIZE;
		real r=rabs(d[0])+rabs(d[1])+rabs(d[2]);
		//fprintf(stderr, "%" RF "g\n", RT(r));
		if(r>=threshold) return 1;
	};
	return 0;
};

int deactivate_slow(const real* grad, int sizep, real threshold) {
	// find maximum
	real max=0;
	for(int n=0; n<SIZE; n++) {		
		if(!ISACTIVE(active,n)) continue;
		for(int p=0; p<sizep; p++) {
			const real* d=grad+3*n+p*3*SIZE;
			real v=rabs(d[0])+rabs(d[1])+rabs(d[2]); 
			if(v>max) max=v;
		};
	};
	max*=threshold;
	//fprintf(stderr, COLOR_BLUE "Max: " COLOR_RESET "%" RF "g\n", RT(max));	
	// turn off slow
	int count=0;
	for(int n=0; n<SIZE; n++) {
		//fprintf(stderr, COLOR_RED "%d: " COLOR_RESET "%d\n", n, count);			
		if(!ISACTIVE(active,n)) continue;
		//fprintf(stderr, COLOR_RED "%d: " COLOR_RESET "%d active\n", n, count);			
		if(is_large(grad+3*n, sizep, max)) count++; 
		else SETPASSIVE(active,n);
	};
	return count;
};

void allocate_nonuniform_field() {
	if(nonuniform_field) return;
	nonuniform_field=ralloc(SIZE*3);
	for(int i=0; i<SIZE; i++) for3(j) 
		nonuniform_field[i*3+j]=magnetic_field[j];
};

void set_tip_field(const real dir[3], const real pos[3]) {
	allocate_nonuniform_field();
	forall(u,x,y,z) {
		real coord[3]; 
		COORDS(u,x,y,z,coord);
		sub3(coord, pos, coord);
		real r=normsq3(coord);
		if(r<1e-8) continue;
		multinv3(r,coord,coord);
		real vec[3];
		cross3(coord, dir, vec);
		int i=INDEX(u,x,y,z);
		for3(j) nonuniform_field[i*3+j]+=vec[j];
	};
};

void prepare_energy_shift(int do_shift) {
	// allocating memory for shifts for every atom in fundamental cell
	zeeman_shift=0;
	energy_shift_per_atom=ralloc(sizeu); // mean energy of S=(0,0,1) over all positions
	anisotropy_shift=0;
	exchange_shift=ralloc(sizeu);
	for(int u=0; u<sizeu; u++) 
		energy_shift_per_atom[u]=exchange_shift[u]=0;
	if(!do_shift) return;
	// zero of energy is set to uniform magnetic fielf with direction ref.
	real ref[3]={0,0,1};
	normalize3(ref);		
	// compute mean zeeman energy of spin=ref
	if(sizen>0 && exchange_constant[0]<0) {

	} else {
		if(nonuniform_field) {
			// TODO: take into account domain/inactive spins
			forall(u,x,y,z) {
				int i=INDEX(u,x,y,z);
				zeeman_shift-=dot3(ref, nonuniform_field+3*i);
			};
			zeeman_shift/=SIZE;
		} else zeeman_shift-=dot3(ref, magnetic_field); 
	};
	// add exchange energy for one atom 
	for(int n=0; n<magnetic_anisotropy_count; n++) {
		real m=dot3(magnetic_anisotropy[n].unit,ref);
		anisotropy_shift-=m*m*magnetic_anisotropy[n].norm;
	};
	for(int n=0;n<sizen;n++) {
		exchange_shift[neighbours[n*5+3]]-=rabs(exchange_constant[n])/2;
		exchange_shift[neighbours[n*5+4]]-=rabs(exchange_constant[n])/2;
	};
	for(int u=0; u<sizeu; u++) 
		energy_shift_per_atom[u]=zeeman_shift+anisotropy_shift+exchange_shift[u];
};

void prepare_dipole_table(real negligible) {
	if(dipole==0) return;
	real sqrt3=rsqrt(3);
	for(int v=0; v<sizeu; v++) {
		real origin[3]; COORDS(v,0,0,0,origin); 
		forall(u,x,y,z) {	
			int lx=x-sizex/2; int ly=y-sizey/2; int lz=z-sizez/2; 
			real dest[3]; COORDS(u,lx,ly,lz,dest);
			sub3(dest,origin,dest);
			real dist=rsqrt(normsq3(dest));
			if(dist==0) continue;
			real alpha=1/(dist*dist*dist);
			if(alpha<negligible) continue;
			if(dipole_count>=dipole_allocated) {
				dipole_allocated=2*dipole_allocated+1;
				dipole_table=(real*)realloc(dipole_table,sizeof(real)*4*dipole_allocated);
				dipole_idx=(int*)realloc(dipole_idx,sizeof(int)*5*dipole_allocated);
			};
			mult3(sqrt3/dist, dest, dipole_table+4*dipole_count+1);
			dipole_table[4*dipole_count]=alpha*dipole;
			dipole_idx[5*dipole_count]=v;
			dipole_idx[5*dipole_count+1]=u;
			dipole_idx[5*dipole_count+2]=lx;
			dipole_idx[5*dipole_count+3]=ly;
			dipole_idx[5*dipole_count+4]=lz;
			dipole_count++;
		};
	};
	fprintf(stderr, "Dipole interaction to be computed on %d neighbours\n", dipole_count);
};
/*
realp skyrmion_minimum_energy() {
	int size=sizeu*sizex*sizey*sizez;
	realp min=0;
	for(int n=0; n<magnetic_anisotropy_count; n++) 
		min-=2*rabs(magnetic_anisotropy[n].norm)*size;
	for(int n=0; n<sizen; n++) {
		min-=2*rsqrt(normsq3(dzyaloshinskii_moriya_vector+3*n))*size;
		min-=2*rabs(exchange_constant[n])*size;
	};
	return min;
};*/

real skyrmion_energy_given_hessian(const real* __restrict__ arg, real* __restrict__ hess_grad) {
	real total_energy=0;
	#pragma omp parallel for collape(4) reduction(+:total_energy)
	forall(u,x,y,z) {
		int i=INDEX(u,x,y,z);
		if(!ISACTIVE(active, i)) continue;
		real* H=nonuniform_field?nonuniform_field+3*i:magnetic_field;
		real G[3]; for3(j) G[j]=hess_grad[3*i+j]/2-H[j];
		for3(j) hess_grad[3*i+j]-=H[j];
		total_energy+=dot3(G, arg+3*i)-energy_shift_per_atom[u];
	};
	return total_energy;
};

void skyrmion_gradient(const real* __restrict__ arg, real* __restrict__ grad, realp* __restrict__ energy) {
	assert(arg); assert(grad);
  	hamiltonian_hessian(arg, grad);
  	if(energy) {
		*energy=skyrmion_energy_given_hessian(arg, grad);
  	} else {
  		subtract_field(grad);  		
  	};
  	
	/*  	
  	if(energy) *energy=-dot(3*SIZE, arg, grad)/2;
  	subtract_field(grad);
  	if(energy) (*energy)+=dot(3*SIZE, arg, grad);
  	*/
};


void projected_gradient(const real* __restrict__ arg, real* __restrict__ grad, realp* __restrict__ energy) {
    skyrmion_gradient(arg, grad, energy);
    project_to_tangent(arg, grad);
};

void hamiltonian_hessian(const real* __restrict__ arg, real* __restrict__ out) {
	// Compute anisotropy part
	real K2[magnetic_anisotropy_count][3]; 
	for(int n=0; n<magnetic_anisotropy_count; n++) for3(j) 
		K2[n][j]=-2*magnetic_anisotropy[n].norm*magnetic_anisotropy[n].unit[j];
	#pragma omp parallel for collapse(4)
	forall(u,x,y,z) {
		int i=INDEX(u,x,y,z);
		for3(j) out[3*i+j]=0.;
		if(!ISACTIVE(active, i)) continue;
		i*=3;
		for(int n=0; n<magnetic_anisotropy_count; n++) {
			real m=dot3(magnetic_anisotropy[n].unit,arg+i);
			for3(j) out[i+j]+=m*K2[n][j];
		};
	};
	// Compute exchange part
	for(int n=0;n<sizen;n++) {
		// local cache
		int s=neighbours[5*n+3], d=neighbours[5*n+4];
		int sx=neighbours[5*n+0], sy=neighbours[5*n+1], sz=neighbours[5*n+2];		
		// Minimum and maximum indices
		int minx, maxx, miny, maxy, minz, maxz; 
		if(boundary_conditions[0]==BC_PERIODIC) { minx=0; maxx=sizex; 
		} else if(sx<0) { minx=-sx; maxx=sizex; 
		} else { maxx=sizex-sx; minx=0; };
		if(boundary_conditions[1]==BC_PERIODIC) { miny=0; maxy=sizey; 
		} else if(sy<0) { miny=-sy; maxy=sizey; 
		} else { maxy=sizey-sy; miny=0; };
		if(boundary_conditions[2]==BC_PERIODIC) { minz=0; maxz=sizez; 
		} else if(sz<0) { minz=-sz; maxz=sizez; 
		} else { maxz=sizez-sz; minz=0; };
		// Compute interaction fo the pair neighbours[n]
		#pragma omp parallel for collapse(3)
		for(int x=minx;x<maxx;x++)for(int y=miny;y<maxy;y++)for(int z=minz;z<maxz;z++) {
			int i1=INDEX(d,(x+sx+sizex)%sizex,(y+sy+sizey)%sizey,(z+sz+sizez)%sizez);
			int i2=INDEX(s,x,y,z);
			if(!ISACTIVE(all_active, i1) || !ISACTIVE(active, i2)) continue;
			i1*=3; i2*=3;
			cross_minus3(dzyaloshinskii_moriya_vector+3*n,arg+i1,out+i2);
			mult_minus3(exchange_constant[n],arg+i1,out+i2);
		};
		#pragma omp parallel for collapse(3)
		for(int x=minx;x<maxx;x++)for(int y=miny;y<maxy;y++)for(int z=minz;z<maxz;z++) {
			int i1=INDEX(d,(x+sx+sizex)%sizex,(y+sy+sizey)%sizey,(z+sz+sizez)%sizez);
			int i2=INDEX(s,x,y,z);
			if(!ISACTIVE(active, i1) || !ISACTIVE(all_active, i2)) continue;
			i1*=3; i2*=3;
			//fprintf(stderr, "%d@ %d %d %d %d -> %d\n",n,d,(x+sx+sizex)%sizex,(y+sy+sizey)%sizey,(z+sz+sizez)%sizez,i1);
			cross_plus3(dzyaloshinskii_moriya_vector+3*n,arg+i2,out+i1);
			mult_minus3(exchange_constant[n],arg+i2,out+i1);
		};
	};
	// Compute dipole interaction
	for(int n=0;n<dipole_count;n++) {
		// local cache
		int s=dipole_idx[5*n+1], d=dipole_idx[5*n+0];
		int sx=dipole_idx[5*n+2], sy=dipole_idx[5*n+3], sz=dipole_idx[5*n+4];
		real* U=dipole_table+4*n+1;
		real alpha=dipole_table[4*n];
		// Minimum and maximum indices
		int minx, maxx, miny, maxy, minz, maxz; 
		if(boundary_conditions[0]==BC_PERIODIC) { minx=0; maxx=sizex; 
		} else if(sx<0) { minx=-sx; maxx=sizex; 
		} else { maxx=sizex-sx; minx=0; };
		if(boundary_conditions[1]==BC_PERIODIC) { miny=0; maxy=sizey; 
		} else if(sy<0) { miny=-sy; maxy=sizey; 
		} else { maxy=sizey-sy; miny=0; };
		if(boundary_conditions[2]==BC_PERIODIC) { minz=0; maxz=sizez; 
		} else if(sz<0) { minz=-sz; maxz=sizez; 
		} else { maxz=sizez-sz; minz=0; };
		// Compute interaction fo the pair neighbours[n]
		#pragma omp parallel for collapse(3)
		for(int x=minx;x<maxx;x++)for(int y=miny;y<maxy;y++)for(int z=minz;z<maxz;z++) {
			int i1=INDEX(d,(x+sx+sizex)%sizex,(y+sy+sizey)%sizey,(z+sz+sizez)%sizez);
			int i2=INDEX(s,x,y,z);
			if(!ISACTIVE(all_active, i1) || !ISACTIVE(active, i2)) continue;
			i1*=3; i2*=3;
			mult_minus3(alpha*dot3(arg+i1,U),U,out+i2);
			mult_minus3(-alpha,arg+i1,out+i2);
		};
		#pragma omp parallel for collapse(3)
		for(int x=minx;x<maxx;x++)for(int y=miny;y<maxy;y++)for(int z=minz;z<maxz;z++) {
			int i1=INDEX(d,(x+sx+sizex)%sizex,(y+sy+sizey)%sizey,(z+sz+sizez)%sizez);
			int i2=INDEX(s,x,y,z);
			if(!ISACTIVE(active, i1) || !ISACTIVE(all_active, i2)) continue;
			i1*=3; i2*=3;
			mult_minus3(alpha*dot3(arg+i2,U),U,out+i1);
			mult_minus3(-alpha,arg+i2,out+i1);
		};
	};	
};

void node_energy(int u, int x, int y, int z, const real* __restrict__ arg, real energy[6]) {
	real anisotropy_energy=0;
	real zeeman_energy=0;
	int i=INDEX(u,x,y,z);
	if(!ISACTIVE(active, i)) {
		for(int k=0; k<6; k++) energy[k]=0.;
		return;
	};
	i*=3;
	if(nonuniform_field) zeeman_energy-=dot3(nonuniform_field+i,arg+i);
	else zeeman_energy-=dot3(magnetic_field,arg+i);
	for(int n=0; n<magnetic_anisotropy_count; n++) {
		real m=dot3(magnetic_anisotropy[n].unit,arg+i);
		anisotropy_energy-=m*m*magnetic_anisotropy[n].norm;
	};
	energy[0]=anisotropy_energy-anisotropy_shift;
	energy[1]=zeeman_energy-zeeman_shift;
	// Compute exchange part
	real dmi_energy=0;
	real exchange_energy=0;
	for(int n=0;n<sizen;n++) {
		// local cache
		int s=neighbours[5*n+3], d=neighbours[5*n+4];
		int sx=neighbours[5*n+0], sy=neighbours[5*n+1], sz=neighbours[5*n+2];		
		// Compute interaction fo the pair neighbours[n]
		if(s==u && (boundary_conditions[0]==BC_PERIODIC || (x+sx<sizex && x+sx>=0))
				&& (boundary_conditions[1]==BC_PERIODIC || (y+sy<sizey && y+sy>=0))
				&& (boundary_conditions[2]==BC_PERIODIC || (z+sz<sizez && z+sz>=0)) ) {
			int i1=INDEX(d,(x+sx+sizex)%sizex,(y+sy+sizey)%sizey,(z+sz+sizez)%sizez);
			int i2=i; // INDEX(s,x,y,z);
			if(!ISACTIVE(all_active, i1)) continue;
			i1*=3; 
			real t[3]; 
			cross3(dzyaloshinskii_moriya_vector+3*n,arg+i1,t);
			dmi_energy-=dot3(t,arg+i2);
			exchange_energy-=exchange_constant[n]*dot3(arg+i1,arg+i2);
		};
		if(d==u && (boundary_conditions[0]==BC_PERIODIC || (x-sx<sizex && x-sx>=0))
				&& (boundary_conditions[1]==BC_PERIODIC || (y-sy<sizey && y-sy>=0))
				&& (boundary_conditions[2]==BC_PERIODIC || (z-sz<sizez && z-sz>=0)) ) {
			int i1=i; // INDEX(d,x,y,z);
			int i2=INDEX(s,(x-sx+sizex)%sizex,(y-sy+sizey)%sizey,(z-sz+sizez)%sizez);
			if(!ISACTIVE(all_active, i2)) continue;
			i2*=3; 
			real t[3]; 
			cross3(dzyaloshinskii_moriya_vector+3*n,arg+i1,t);
			dmi_energy-=dot3(t,arg+i2);
			exchange_energy-=exchange_constant[n]*dot3(arg+i1,arg+i2);
		};		
	};
	energy[3]=dmi_energy;
	energy[2]=exchange_energy-exchange_shift[u];
	// Compute dipole interaction
	real dipole_energy=0;	
	for(int n=0;n<dipole_count;n++) {
		// local cache
		int s=dipole_idx[5*n+1], d=dipole_idx[5*n+0];
		int sx=dipole_idx[5*n+2], sy=dipole_idx[5*n+3], sz=dipole_idx[5*n+4];
		real* U=dipole_table+4*n+1;
		real alpha=dipole_table[4*n];
		
		// Compute interaction fo the pair neighbours[n]
		if(s==u && (boundary_conditions[0]==BC_PERIODIC || (x+sx<sizex && x+sx>=0))
				&& (boundary_conditions[1]==BC_PERIODIC || (y+sy<sizey && y+sy>=0))
				&& (boundary_conditions[2]==BC_PERIODIC || (z+sz<sizez && z+sz>=0)) ) {
			int i1=INDEX(d,(x+sx+sizex)%sizex,(y+sy+sizey)%sizey,(z+sz+sizez)%sizez);
			int i2=i;
			if(!ISACTIVE(all_active, i1)) continue;
			i1*=3; 
			dipole_energy-=alpha*(dot3(arg+i1,U)*dot3(arg+i2,U)-dot3(arg+i2,arg+i1));
		};
		if(d==u && (boundary_conditions[0]==BC_PERIODIC || (x-sx<sizex && x-sx>=0))
				&& (boundary_conditions[1]==BC_PERIODIC || (y-sy<sizey && y-sy>=0))
				&& (boundary_conditions[2]==BC_PERIODIC || (z-sz<sizez && z-sz>=0)) ) {
			int i1=i; //INDEX(d,x,y,z);
			int i2=INDEX(s,(x-sx+sizex)%sizex,(y-sy+sizey)%sizey,(z-sz+sizez)%sizez);
			if(!ISACTIVE(all_active, i2)) continue;
			i2*=3;
			dipole_energy-=alpha*(dot3(arg+i1,U)*dot3(arg+i2,U)-dot3(arg+i2,arg+i1));
		};
	};	
	energy[4]=dipole_energy;
	energy[5]=energy[0]+energy[1]+energy[2]+energy[3]+energy[4];
};

// energy[0] - anisotropy part
// energy[1] - zeeman part (external field)
// energy[2] - exchange part
// energy[3] - D-M energy
// energy[4] - dipole energy
// energy[5] - total energy
void skyrmion_energy(const real* __restrict__ arg, realp energy[6]) {
	//for(int j=0;j<5;j++) energy[j]=0;
	// Compute anisotropy part
	energy[0]=0;
	energy[1]=0;
	int countu[sizeu]; for(int u=0; u<sizeu; u++) countu[u]=0;
	// TODO: fix reduction on vector countu 
	#pragma omp parallel for collapse(4) reduction(+:anisotropy_energy,zeeman_energy,countu)
	forall(u,x,y,z) {
		realp anisotropy_energy=0;
		realp zeeman_energy=0;
		int i=INDEX(u,x,y,z);
		if(!ISACTIVE(active, i)) continue;
		countu[u]++;
		i*=3;
		for(int n=0; n<magnetic_anisotropy_count; n++) {
			real m=dot3(magnetic_anisotropy[n].unit,arg+i);
			anisotropy_energy-=m*m*magnetic_anisotropy[n].norm;
		};
		if(nonuniform_field) zeeman_energy-=dot3(nonuniform_field+i,arg+i);
		else zeeman_energy-=dot3(magnetic_field,arg+i);
		energy[0]+=anisotropy_energy-anisotropy_shift;
		energy[1]+=zeeman_energy-zeeman_shift;
	};
	int count=0;
	for(int u=0; u<sizeu; u++) count+=countu[u];

	// Compute exchange part
	realp dmi_energy=0;
	realp exchange_energy[sizeu];
	for(int u=0; u<sizeu; u++) { exchange_energy[u]=0; };
	for(int n=0;n<sizen;n++) {
		// local cache
		int s=neighbours[5*n+3], d=neighbours[5*n+4];
		int sx=neighbours[5*n+0], sy=neighbours[5*n+1], sz=neighbours[5*n+2];		
		// Minimum and maximum indices
		int minx, maxx, miny, maxy, minz, maxz; 
		if(boundary_conditions[0]==BC_PERIODIC) { minx=0; maxx=sizex; 
		} else if(sx<0) { minx=-sx; maxx=sizex; 
		} else { maxx=sizex-sx; minx=0; };
		if(boundary_conditions[1]==BC_PERIODIC) { miny=0; maxy=sizey; 
		} else if(sy<0) { miny=-sy; maxy=sizey; 
		} else { maxy=sizey-sy; miny=0; };
		if(boundary_conditions[2]==BC_PERIODIC) { minz=0; maxz=sizez; 
		} else if(sz<0) { minz=-sz; maxz=sizez; 
		} else { maxz=sizez-sz; minz=0; };
		// Compute interaction fo the pair neighbours[n]
		#pragma omp parallel for collapse(3) reduction(+:dmi_energy,exchange_energy)
		for(int x=minx;x<maxx;x++)for(int y=miny;y<maxy;y++)for(int z=minz;z<maxz;z++) {
			int i1=INDEX(d,(x+sx+sizex)%sizex,(y+sy+sizey)%sizey,(z+sz+sizez)%sizez);
			int i2=INDEX(s,x,y,z);
			if(!(  (ISACTIVE(active, i1) && ISACTIVE(all_active, i2))
				|| (ISACTIVE(all_active, i1) && ISACTIVE(active, i2))
				)
			  ) continue;
			i1*=3; i2*=3;
			real t[3]; 
			cross3(dzyaloshinskii_moriya_vector+3*n,arg+i1,t);
			dmi_energy-=dot3(t,arg+i2);
			realp exen=exchange_constant[n]*dot3(arg+i1,arg+i2)/2;
			exchange_energy[s]-=exen; 
			exchange_energy[d]-=exen; 
		};
	};
	energy[3]=dmi_energy;
	energy[2]=0;
	for(int u=0; u<sizeu; u++) energy[2]+=exchange_energy[u]-exchange_shift[u]*countu[u];
	// Compute dipole interaction
	realp dipole_energy=0;	
	for(int n=0;n<dipole_count;n++) {
		// local cache
		int s=dipole_idx[5*n+1], d=dipole_idx[5*n+0];
		int sx=dipole_idx[5*n+2], sy=dipole_idx[5*n+3], sz=dipole_idx[5*n+4];
		real* U=dipole_table+4*n+1;
		real alpha=dipole_table[4*n];
		// Minimum and maximum indices
		int minx, maxx, miny, maxy, minz, maxz; 
		if(boundary_conditions[0]==BC_PERIODIC) { minx=0; maxx=sizex; 
		} else if(sx<0) { minx=-sx; maxx=sizex; 
		} else { maxx=sizex-sx; minx=0; };
		if(boundary_conditions[1]==BC_PERIODIC) { miny=0; maxy=sizey; 
		} else if(sy<0) { miny=-sy; maxy=sizey; 
		} else { maxy=sizey-sy; miny=0; };
		if(boundary_conditions[2]==BC_PERIODIC) { minz=0; maxz=sizez; 
		} else if(sz<0) { minz=-sz; maxz=sizez; 
		} else { maxz=sizez-sz; minz=0; };
		// Compute interaction fo the pair neighbours[n]
		#pragma omp parallel for collapse(3) reduction(+:dipole_energy)
		for(int x=minx;x<maxx;x++)for(int y=miny;y<maxy;y++)for(int z=minz;z<maxz;z++) {
			int i1=INDEX(d,(x+sx+sizex)%sizex,(y+sy+sizey)%sizey,(z+sz+sizez)%sizez);
			int i2=INDEX(s,x,y,z);
			if(!(  (ISACTIVE(active, i1) && ISACTIVE(all_active, i2))
				|| (ISACTIVE(all_active, i1) && ISACTIVE(active, i2))
				)
			  ) continue;
			i1*=3; i2*=3;
			dipole_energy-=alpha*(dot3(arg+i1,U)*dot3(arg+i2,U)-dot3(arg+i2,arg+i1));
		};
	};	
	energy[4]=dipole_energy;	
	energy[5]=energy[0]+energy[1]+energy[2]+energy[3]+energy[4];
};

void subtract_field(real* __restrict__ inout) {
	if(nonuniform_field) {
		#pragma omp parallel for collapse(4)
		forall(u,x,y,z) {
			int i=INDEX(u,x,y,z);
			if(!ISACTIVE(active, i)) continue;
			i*=3;
			for3(j) inout[i+j]-=nonuniform_field[i+j];
			};
	} else {
		#pragma omp parallel for collapse(4)
		forall(u,x,y,z) for3(j) {
			int i=INDEX(u,x,y,z);
			if(!ISACTIVE(active, i)) continue;
			inout[i*3+j]-=magnetic_field[j];
		}
	};
};

void set_to_field(real* __restrict__ out) {
	real field[3]={0,0,1};
	#pragma omp parallel for collapse(4)
	forall(u,x,y,z) {
		int i=INDEX(u,x,y,z);
		if(!ISACTIVE(active, i)) for3(j) out[i*3+j]=0;
		else for3(j) out[i*3+j]=field[j];
	};
};

// Normalize vector field so every vector has unit length 
real normalize(real* __restrict__ a) {
	real sum=0;
	#pragma omp parallel for collapse(4) reduction(+:sum)
	forall(u,x,y,z) {
		int i=INDEX(u,x,y,z);
		if(!ISACTIVE(active, i)) continue;
		sum+=normalize3(a+i*3);
	};
	return sum;
};

realp seminormalize(real factor, real* __restrict__ a) {
	realp sum=0;
	#pragma omp parallel for collapse(4) reduction(+:sum)
	forall(u,x,y,z) {
		int i=INDEX(u,x,y,z);
		if(!ISACTIVE(active, i)) continue;
		sum+=seminormalize3(factor,a+i*3);
	};
	return sum;
};

// Project vector field 't' to tangent space of unit length vector field 'a'
void project_to_tangent(const real* __restrict__ a, real* __restrict__ b) {
	#pragma omp parallel for collapse(4)
	forall(u,x,y,z) {
		int i=INDEX(u,x,y,z);
		if(!ISACTIVE(active, i)) { for3(j) b[3*i+j]=0; continue; };
		i*=3; tangent3(a+i,b+i);
	};
};

// C:x->(<x|P_j x>/2-1/2)_j
void skyrmion_constrain(const real* __restrict__ a, real* __restrict__ r) {
	#pragma omp parallel for collapse(4)	
  	forall(u,x,y,z) {
  		int i=INDEX(u,x,y,z);
  		r[i]=ISACTIVE(active, i)?(normsq3(a+3*i)-1.)/2.:0.;
  	};
};

// D:x,u,r->r+sum_l u_l P_j x
void skyrmion_constrain_gradient(const real* __restrict__ a, const real* __restrict__ lambda, real* __restrict__ r) {
	#pragma omp parallel for collapse(4)
  	forall(u,x,y,z) {
	  	int i=INDEX(u,x,y,z);
	  	if(!ISACTIVE(active, i)) continue;
	  	mult_plus3(lambda[i],a+3*i,r+3*i);
  	};
};

// P:x,y->(<x|P_j y>)_l
void skyrmion_constrain_adjucent(const real* __restrict__ a, const real* __restrict__ b, real* __restrict__ r) {
	#pragma omp parallel for collapse(4)	
  	forall(u,x,y,z) {
	  	int i=INDEX(u,x,y,z);
	  	r[i]=ISACTIVE(active, i)?dot3(a+3*i,b+3*i):0.;
  	};
};

void skyrmion_middle(const real* __restrict__ a, const real* __restrict__ b, real* __restrict__ r) {
	#pragma omp parallel for collapse(4)	
  	forall(u,x,y,z) {
	  	int i=INDEX(u,x,y,z);
	  	if(!ISACTIVE(active, i))  { for3(j) r[3*i+j]=0; continue; };
	  	middle3(a+3*i,b+3*i,r+3*i);
  	};
};

// calculate point in between of b and c by interpolation of curve given by
// points a,b,c and d in that order
void skyrmion_middle_fourth_order(const real* __restrict__ a, const real* __restrict__ b, const real* __restrict__ c, const real* __restrict__ d, real* __restrict__ r) {
	#pragma omp parallel for collapse(4)	
  	forall(u,x,y,z) {
	  	int i=INDEX(u,x,y,z);
	  	if(!ISACTIVE(active, i))  { for3(j) r[3*i+j]=0; continue; };
	  	middle_fourth_order3(a+3*i,b+3*i,c+3*i,d+3*i,r+3*i);
  	};
};

void skyrmion_middle_third_order(const real* __restrict__ a, const real* __restrict__ b, const real* __restrict__ c, real* __restrict__ r) {
	#pragma omp parallel for collapse(4)	
  	forall(u,x,y,z) {
	  	int i=INDEX(u,x,y,z);
	  	if(!ISACTIVE(active, i))  { for3(j) r[3*i+j]=0; continue; };
	  	middle_third_order3(a+3*i,b+3*i,c+3*i,r+3*i);
  	};
};
	
void skyrmion_geodesic_rec(real noise, real* p, int n, int m) {
	int size=SIZE*3;
	if(n+1>=m) return;
	int k=(n+m)/2; // Точка разбиения
	// Находим середину на прямой между n и m и проецируем на сферы
	skyrmion_middle(p+n*size,p+m*size,p+k*size);
	if(noise!=0) {
		add_random_vector(noise*(m-n),size,p+k*size,p+k*size);
		normalize(p+k*size);
	};
	// Заполняем пробелы
	skyrmion_geodesic_rec(noise,p,n,k);
	skyrmion_geodesic_rec(noise,p,k,m);
}

void skyrmion_geodesic(real noise, int sizep, real* p) { 
	skyrmion_geodesic_rec(noise, p, 0, sizep-1); 
}

void two_point_tangent0(const real* __restrict__ a, const real* __restrict__ b, real* __restrict__ r) {
	#pragma omp parallel for collapse(4)	
	forall(u,x,y,z) {	
		int i=INDEX(u,x,y,z);
		if(!ISACTIVE(active, i))  { for3(j) r[3*i+j]=0; continue; };
		i*=3;
		sub3(b+i,a+i,r+i); 
		tangent3(a+i,r+i); 
		//normalize3(r+i);
	};
};

void two_point_tangent1(const real* __restrict__ a, const real* __restrict__ b, real* __restrict__ r) {
	#pragma omp parallel for collapse(4)	
	forall(u,x,y,z) {	
		int i=INDEX(u,x,y,z);
		if(!ISACTIVE(active, i))  { for3(j) r[3*i+j]=0; continue; };
		i*=3;
		sub3(b+i,a+i,r+i); 
		tangent3(b+i,r+i); 
		//normalize3(r+i);
	};
};

// tangent r to path defined by three consequative points a,b,c
void three_point_tangent(const real* __restrict__ a, const real* __restrict__ b, const real* __restrict__ c, real* __restrict__ r) {
	#pragma omp parallel for collapse(4)	
	forall(u,x,y,z) {	
		int i=INDEX(u,x,y,z);
		if(!ISACTIVE(active, i))  { for3(j) r[3*i+j]=0; continue; };
		i*=3;
		sub3(c+i,a+i,r+i); 
		tangent3(b+i,r+i); 
		//normalize3(r+i);
	};
};

void three_point_tangent_stable(real ea, real eb, real ec, const real* __restrict__ a, const real* __restrict__ b, const real* __restrict__ c, real* __restrict__ r) {
	if(ea<eb && eb<ec) {
		#pragma omp parallel for collapse(4)	
		forall(u,x,y,z) {	
			int i=INDEX(u,x,y,z);
			if(!ISACTIVE(active, i))  { for3(j) r[3*i+j]=0; continue; };
			i*=3;
			sub3(c+i,b+i,r+i); 
			tangent3(b+i,r+i); 
		};
	} else if(ea>eb && eb>ec) {
		#pragma omp parallel for collapse(4)	
		forall(u,x,y,z) {	
			int i=INDEX(u,x,y,z);
			if(!ISACTIVE(active, i))  { for3(j) r[3*i+j]=0; continue; };
			i*=3;
			sub3(b+i,a+i,r+i); 
			tangent3(b+i,r+i); 
		};
	} else {
		real w1=rabs(eb-ea); real w2=rabs(ec-eb);
		if(w1>w2) { real t=w1; w1=w2; w2=t; }; // w1=min(|eb-ea|,|ec-eb|), w2=max(...)
		if(ec<ea) { real t=w1; w1=w2; w2=t; }; // w1=min if ec>ea  w1=max otherwise
		#pragma omp parallel for collapse(4)	
		forall(u,x,y,z) {	
			int i=INDEX(u,x,y,z);
			if(!ISACTIVE(active, i))  { for3(j) r[3*i+j]=0; continue; };
			i*=3;
			real t1[3],t2[3];
			sub3(b+i,a+i,t1); 
			sub3(c+i,b+i,t2);
			for3(j) (r+i)[j]=w1*t1[j]+w2*t2[j];
			tangent3(b+i,r+i); 
		};
	}
};

// tangent r to path defined by three consequative points a,b,c
void three_point_tangent_mean(const real* __restrict__ a, const real* __restrict__ b, const real* __restrict__ c, real* __restrict__ r) {
	#pragma omp parallel for collapse(4)	
	forall(u,x,y,z) {	
		real t1[3],t2[3];
		int i=INDEX(u,x,y,z);
		if(!ISACTIVE(active, i))  { for3(j) r[3*i+j]=0; continue; };
		i*=3;
		sub3(b+i,a+i,t1); 
		sub3(c+i,b+i,t2); 
		real l1=1./rsqrt(normsq3(t1)); real l2=1./rsqrt(normsq3(t2));
		for3(j) (r+i)[j]=t1[j]*l1+t2[j]*l2;
		tangent3(b+i,r+i); 
		//normalize3(r+i);
	};
};


/*
// r moved along b-a to satisfy |r-a|=|r-b|
void three_point_equalize(const real* __restrict__ a, const real* __restrict__ b, real* __restrict__ r) {
	#pragma omp parallel for collapse(4)	
	forall(u,x,y,z) {	
		real vec[3];
		real ab,ar,br;
		int i=INDEX(u,x,y,z)*3;
		sub3(r+i,a+i,vec); ar=normsq3(vec);
		sub3(r+i,b+i,vec); br=normsq3(vec);
		sub3(b+i,a+i,vec); ab=normsq3(vec);
		mult_plus3((br-ar)/ab*0.5,vec,r+i);
		normalize3(r+i);
	};
};

void three_point_equalizer(const real* __restrict__ a, const real* __restrict__ c, const real* __restrict__ b, real* __restrict__ r) {
	#pragma omp parallel for collapse(4)	
	forall(u,x,y,z) {	
		real vec[3];
		real ab,ar,br;
		int i=INDEX(u,x,y,z)*3;
		sub3(c+i,a+i,vec); ar=normsq3(vec);
		sub3(c+i,b+i,vec); br=normsq3(vec);
		sub3(b+i,a+i,vec); ab=normsq3(vec);
		copy3(c+i,r+i);
		if(ab>10*EPSILON) mult_plus3((br-ar)/ab*0.5,vec,r+i);
		normalize3(r+i);
	};
};
*/
void append_skyrmion(const real center[3], real distance, 
	real winding_rho, real winding_phi, real winding_z, 
	real rotation_rho, real rotation_phi, real rotation_z, 
	real z_rot_rho, real z_rot_phi, real z_rot_z, 
	real* __restrict__ data)
{
	real field[3]={0,0,1}; 
	#pragma omp parallel for collapse(4)	
	forall(u,x,y,z) {	
		int i=INDEX(u,x,y,z);
		if(!ISACTIVE(active, i)) continue;
		i*=3;
		real vec[3]; COORDS(u,x,y,z,vec);
		sub3(vec,center,vec);
		real hei=vec[2];
		real distsq=vec[0]*vec[0]+vec[1]*vec[1];
		if(distsq>distance*distance-hei*hei) continue;
		if(distsq==0) {
			mult3(-1,data+i,data+i);
			continue;
		};
		vec[2]=0; normalize3(vec);
		real phi=ratan2(vec[1],vec[0]);
		// First rotation is around vec
		real dist=rsqrt(distsq)/distance; hei/=distance; 
		dist/=rsqrt(1-hei*hei); dist=1-dist; 
		dist*=M_PI_2; hei*=M_PI_2;
		real sinalpha, cosalpha;
		rsincos(dist*winding_rho+phi*winding_phi+hei*winding_z,&sinalpha,&cosalpha); 
		real q1[4]; q1[0]=cosalpha; mult3(sinalpha,vec,q1+1);
		// Second rotation is in the plane containing vec and field
		rsincos(dist*rotation_rho+phi*rotation_phi+hei*rotation_z,&sinalpha,&cosalpha); 
		real q2[4]; q2[0]=cosalpha; cross3(field,vec,q2+1);
		mult3(sinalpha,q2+1,q2+1);
		// Third rotation is around z axis
		rsincos(dist*z_rot_rho+phi*z_rot_phi+hei*z_rot_z,&sinalpha,&cosalpha); 
		real q3[4]={cosalpha,0,0,sinalpha}; 
		// combined rotation
		real qt[4]; quaternion_product(q2,q1,qt);
		real q[4]; quaternion_product(q3,qt,q);
		// applying rotation
		q1[0]=0; copy3(data+i, q1+1);
		quaternion_product(q,q1,q2);
		q[1]=-q[1]; q[2]=-q[2]; q[3]=-q[3];
		quaternion_product(q2,q,q1);
		//fprintf(stderr,"Error %" RF "g\n",rsqrt(normsq3(q1+1))-1);
		//assert(rabs(rsqrt(normsq3(q1+1))-1)<1e-8);
		copy3(q1+1,data+i);
	};
}


// axis: 0=x, 1=y, 2=z
void group_generator(const real* __restrict__ spins, 
int axis, real* __restrict__ gen) {
	int N;
	switch(axis) {
		case 0: N=sizex; break;
		case 1: N=sizey; break;
		case 2: N=sizez; break;
		default: 
			fprintf(stderr, "Wrong axis: %d not in 0..2\n",axis);
			exit(1);
	};
// checking parameters
	if(N%2==0) {
		fprintf(stderr, COLOR_BOLD COLOR_RED "Error:" COLOR_RESET 
			"Axis %d has even length %d\n",axis,N);
		exit(1);
	};
// computing kernel	
	real* kernel=ralloc(N);
	real mult=M_PI/N; 
	kernel[0]=0;
	for(int n=1; n<N; n++) {
		real s,c; rsincos(mult*n,&s,&c);
		kernel[n]=mult/s;
		if(n%2==1) kernel[n]=-kernel[n];
	};
// doing covolution
	forall(u,x,y,z) for3(j) gen[j+3*INDEX(u,x,y,z)]=0;
	for(int n=0; n<N; n++) {
		real a=kernel[n]; 
		if(a==0) continue;
		forall(u,x,y,z) {
			int i=3*INDEX(u,x,y,z); int k;
			switch(axis) {
				case 0: k=3*INDEX(u,n>=x?n-x:n-x+N,y,z); break;
				case 1: k=3*INDEX(u,x,n>=y?n-y:n-y+N,z); break;
				case 2: k=3*INDEX(u,x,y,n>=z?n-z:n-z+N); break;
				default: assert(0);
			};
			for3(j) gen[j+i]+=a*spins[j+k];
		};
	};

// finilizing
	free(kernel);
};

void skyrmion_afm(real* __restrict__ a) {
	forall(u,x,y,z) {
		int n=INDEX(u,x,y,z); int i=3*n;
		if(!ISACTIVE(active, n)) {
			for3(j) a[i+j]=0;
		} else if(((x+y+z)*sizeu+u)%2) { 
			for3(j) a[i+j]*=-1;
		};
	}
};

void skyrmion_random(real* __restrict__ a) {
	forall(u,x,y,z) {
		int i=INDEX(u,x,y,z);
		if(ISACTIVE(active, i)) {
			i*=3;
			for3(j) a[i+j]=random_real()*2-1;
		} else {
			i*=3;
			for3(j) a[i+j]=0;
		};
	}
	normalize(a);
};
