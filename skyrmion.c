#define _GNU_SOURCE
#include <math.h>
#include <stdlib.h>

#include "vector.h"
#include "skyrmion.h"

// Physical parameters
real magnetic_field[3]={NAN,NAN,NAN};
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
real magnetic_anisotropy_norm=NAN;
real magnetic_anisotropy_unit[3]={NAN,NAN,NAN};
// Exchange constant J
real* exchange_constant=NULL;
// Dzyaloshinskii Moriya vector for every pair of atoms
real* dzyaloshinskii_moriya_vector=NULL;

void hamiltonian_hessian(const real* restrict arg, real* restrict out) {
	// Compute anisotropy part
	real K2[3]; for3(j) K2[j]=-2*magnetic_anisotropy_norm*magnetic_anisotropy_unit[j];
	forall(u,x,y,z) {
		real m=dot3(magnetic_anisotropy_unit,arg+3*INDEX(u,x,y,z));
		for3(j) out[3*INDEX(u,x,y,z)+j]=m*K2[j];
	};
	// Comput exchange part
	for(int n=0;n<sizen;n++) {
		// local cache
		int d=neighbours[5*n+3], s=neighbours[5*n+4];
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
		for(int x=minx;x<maxx;x++)for(int y=miny;y<maxy;y++)for(int z=minz;z<maxz;z++) {
			cross_minus3(dzyaloshinskii_moriya_vector+3*n,arg+INDEX(s,(x+sx)%sizex,(y+sy)%sizey,(z+sz)%sizez)*3,out+INDEX(d,x,y,z)*3);
			mult_minus3(exchange_constant[n],arg+INDEX(s,(x+sx)%sizex,(y+sy)%sizey,(z+sz)%sizez)*3,out+INDEX(d,x,y,z)*3);
		};
		for(int x=minx;x<maxx;x++)for(int y=miny;y<maxy;y++)for(int z=minz;z<maxz;z++) {
			cross_plus3(dzyaloshinskii_moriya_vector+3*n,arg+INDEX(s,x,y,z)*3,out+INDEX(d,(x+sx)%sizex,(y+sy)%sizey,(z+sz)%sizez)*3);
			mult_minus3(exchange_constant[n],arg+INDEX(s,x,y,z)*3,out+INDEX(d,(x+sx)%sizex,(y+sy)%sizey,(z+sz)%sizez)*3);
		};
	};
};

void subtract_field(real* restrict inout) {
	forall(u,x,y,z) for3(j) inout[INDEX(u,x,y,z)*3+j]-=magnetic_field[j];
};

// Normalize vector field so every vector has unit length 
void normalize(real* restrict a) {
	forall(u,x,y,z) {
		normalize3(a+INDEX(u,x,y,z)*3);
	};
};

// Project vector field 't' to tangent space of unit length vector field 'a'
void project_to_tangent(const real* restrict a, real* restrict b) {
	forall(u,x,y,z) {
		tangent3(a+INDEX(u,x,y,z)*3,b+INDEX(u,x,y,z)*3);
	};
};

// C:x->(<x|P_j x>/2-1/2)_j
void skyrmion_constrain(const real* restrict a, real* restrict r) {
  forall(u,x,y,z) {
  	int i=INDEX(u,x,y,z);
  	r[i]=(normsq3(a+3*i)-1.)/2.;
  };
};

// D:x,u,r->r+sum_l u_l P_j x
void skyrmion_constrain_gradient(const real* restrict a, const real* restrict lambda, real* restrict r) {
  forall(u,x,y,z) {
  	int i=INDEX(u,x,y,z);
  	mult_plus3(lambda[i],a+3*i,r+3*i);
  };
};

// P:x,y->(<x|P_j y>)_l
void skyrmion_constrain_adjucent(const real* restrict a, const real* restrict b, real* restrict r) {
  forall(u,x,y,z) {
  	int i=INDEX(u,x,y,z);
  	r[i]=dot3(a+3*i,b+3*i);
  };
};







void fourier_table(const real* restrict angles, real* restrict table) {
	forall(u,x,y,z) {
		int i=INDEX(u,x,y,z);
		rsincos(angles[2*i+0],table+4*i+0,table+4*i+1);
		rsincos(angles[2*i+1],table+4*i+2,table+4*i+3);
	};
};

// transform angles to vector on sphere
// (x,y,z)=R(phi,theta)
void angles_to_vector(const real* restrict table, const real* restrict angles, real* restrict vectors) {
	forall(u,x,y,z) {
		int i=INDEX(u,x,y,z);
		real sphi=table[4*i+0];
		real cphi=table[4*i+1];
		real stheta=table[4*i+2];
		real ctheta=table[4*i+3];
		vectors[3*i+0]=stheta*cphi;
		vectors[3*i+1]=stheta*sphi;
		vectors[3*i+2]=ctheta;
	};
};

// project vectors to tangent space on sphere
// (x,y,z)*grad_{x,y,z} R(phi,theta)
/*void tangent_vector_to_angles(const field table[4], const field vector[3], field angles[2]) {
	forall(u,x,y,z) {
		real sphi=table[u][x][y][z][0];
		real cphi=table[u][x][y][z][1];
		real stheta=table[u][x][y][z][2];
		real ctheta=table[u][x][y][z][3];
		angles[u][x][y][z][0]=;
		angles[u][x][y][z][1]=;
	};	
};*/	