#define _GNU_SOURCE
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "vector.h"
#include "skyrmion.h"
#include "debug.h"

// Physical parameters
real magnetic_field[3]={NAN,NAN,NAN};
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
real magnetic_anisotropy_norm=NAN;
real magnetic_anisotropy_unit[3]={NAN,NAN,NAN};
// Exchange constant J
real* exchange_constant=NULL;
// dipole interaction constant
real dipole=0;
// Dzyaloshinskii Moriya vector for every pair of atoms
real* dzyaloshinskii_moriya_vector=NULL;
real* initial_state=NULL;
int initial_states_count=0;

int dipole_count=0;
int dipole_allocated=0;
int* dipole_idx=NULL; // <source atom> <dest atom> <dest x> <dest y> <dest z>
real* dipole_table=NULL; // <multiplier> <vector x> <vector y> <vector z>


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

real skyrmion_minimum_energy() {
	int size=sizeu*sizex*sizey*sizez;
	real min=-2*rabs(magnetic_anisotropy_norm)*size;
	for(int n=0; n<sizen;n++) {
		min-=2*rsqrt(normsq3(dzyaloshinskii_moriya_vector+3*n))*size;
		min-=2*rabs(exchange_constant[n])*size;
	};
	return min;
};

void hamiltonian_hessian(const real* restrict arg, real* restrict out) {
	// Compute anisotropy part
	real K2[3]; for3(j) K2[j]=-2*magnetic_anisotropy_norm*magnetic_anisotropy_unit[j];
	#pragma omp parallel for collapse(4)
	forall(u,x,y,z) {
		int i=3*INDEX(u,x,y,z);
		real m=dot3(magnetic_anisotropy_unit,arg+i);
		for3(j) out[i+j]=m*K2[j];
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
			int i1=INDEX(d,(x+sx+sizex)%sizex,(y+sy+sizey)%sizey,(z+sz+sizez)%sizez)*3;
			int i2=INDEX(s,x,y,z)*3;
			cross_minus3(dzyaloshinskii_moriya_vector+3*n,arg+i1,out+i2);
			mult_minus3(exchange_constant[n],arg+i1,out+i2);
		};
		#pragma omp parallel for collapse(3)
		for(int x=minx;x<maxx;x++)for(int y=miny;y<maxy;y++)for(int z=minz;z<maxz;z++) {
			int i1=INDEX(d,(x+sx+sizex)%sizex,(y+sy+sizey)%sizey,(z+sz+sizez)%sizez)*3;
			int i2=INDEX(s,x,y,z)*3;
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
			int i1=INDEX(d,(x+sx+sizex)%sizex,(y+sy+sizey)%sizey,(z+sz+sizez)%sizez)*3;
			int i2=INDEX(s,x,y,z)*3;
			mult_minus3(alpha*dot3(arg+i1,U),U,out+i2);
			mult_minus3(-alpha,arg+i1,out+i2);
		};
		#pragma omp parallel for collapse(3)
		for(int x=minx;x<maxx;x++)for(int y=miny;y<maxy;y++)for(int z=minz;z<maxz;z++) {
			int i1=INDEX(d,(x+sx+sizex)%sizex,(y+sy+sizey)%sizey,(z+sz+sizez)%sizez)*3;
			int i2=INDEX(s,x,y,z)*3;
			mult_minus3(alpha*dot3(arg+i2,U),U,out+i1);
			mult_minus3(-alpha,arg+i2,out+i1);
		};
	};	
};

// energy[0] - anisotropy part
// energy[1] - zeeman part (external field)
// energy[2] - exchange part
// energy[3] - D-M energy
// energy[4] - dipole energy
void skyrmion_energy(const real* restrict arg, real energy[5]) {
	//for(int j=0;j<5;j++) energy[j]=0;
	// Compute anisotropy part
	real anisotropy_energy=0;
	real zeeman_energy=0;
	#pragma omp parallel for collapse(3) reduction(+:anisotropy_energy,zeeman_energy)
	forall(u,x,y,z) {
		int i=3*INDEX(u,x,y,z);
		real m=dot3(magnetic_anisotropy_unit,arg+i);
		anisotropy_energy-=m*m;
		if(nonuniform_field) zeeman_energy-=dot3(nonuniform_field+i,arg+i);
		else zeeman_energy-=dot3(magnetic_field,arg+i);
	};
	energy[0]=anisotropy_energy*magnetic_anisotropy_norm;
	energy[1]=zeeman_energy;
	// Compute exchange part
	real dmi_energy=0;
	real exchange_energy=0;
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
			int i1=INDEX(d,(x+sx+sizex)%sizex,(y+sy+sizey)%sizey,(z+sz+sizez)%sizez)*3;
			int i2=INDEX(s,x,y,z)*3;
			real t[3]; 
			cross3(dzyaloshinskii_moriya_vector+3*n,arg+i1,t);
			dmi_energy-=dot3(t,arg+i2);
			exchange_energy-=exchange_constant[n]*dot3(arg+i1,arg+i2);
		};
	};
	energy[3]=dmi_energy;
	energy[2]=exchange_energy;
	// Compute dipole interaction
	real dipole_energy=0;	
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
			int i1=INDEX(d,(x+sx+sizex)%sizex,(y+sy+sizey)%sizey,(z+sz+sizez)%sizez)*3;
			int i2=INDEX(s,x,y,z)*3;
			dipole_energy-=alpha*(dot3(arg+i1,U)*dot3(arg+i2,U)-dot3(arg+i2,arg+i1));
		};
		energy[4]=dipole_energy;
	};		
};

void subtract_field(real* restrict inout) {
	if(nonuniform_field) {
		#pragma omp parallel for collapse(4)
		forall(u,x,y,z) {
			int i=INDEX(u,x,y,z)*3;
			for3(j) inout[i+j]-=nonuniform_field[i+j];
			};
	} else {
		#pragma omp parallel for collapse(4)
		forall(u,x,y,z) for3(j) 
			inout[INDEX(u,x,y,z)*3+j]-=magnetic_field[j];
	};
};

void set_to_field(real* restrict out) {
	real field[3]={0,0,1};
	#pragma omp parallel for collapse(4)
	forall(u,x,y,z) for3(j) out[INDEX(u,x,y,z)*3+j]=field[j];
};


// Normalize vector field so every vector has unit length 
void normalize(real* restrict a) {
	#pragma omp parallel for collapse(4)
	forall(u,x,y,z) normalize3(a+INDEX(u,x,y,z)*3);
};

real seminormalize(real factor, real* restrict a) {
	real sum=0;
	#pragma omp parallel for collapse(4) reduction(+:sum)
	forall(u,x,y,z) sum+=seminormalize3(factor,a+INDEX(u,x,y,z)*3);
	return sum;
};

// Project vector field 't' to tangent space of unit length vector field 'a'
void project_to_tangent(const real* restrict a, real* restrict b) {
	#pragma omp parallel for collapse(4)
	forall(u,x,y,z) {
		int i=INDEX(u,x,y,z)*3;
		tangent3(a+i,b+i);
	};
};

// C:x->(<x|P_j x>/2-1/2)_j
void skyrmion_constrain(const real* restrict a, real* restrict r) {
	#pragma omp parallel for collapse(4)	
  	forall(u,x,y,z) {
  		int i=INDEX(u,x,y,z);
  		r[i]=(normsq3(a+3*i)-1.)/2.;
  	};
};

// D:x,u,r->r+sum_l u_l P_j x
void skyrmion_constrain_gradient(const real* restrict a, const real* restrict lambda, real* restrict r) {
	#pragma omp parallel for collapse(4)
  	forall(u,x,y,z) {
	  	int i=INDEX(u,x,y,z);
	  	mult_plus3(lambda[i],a+3*i,r+3*i);
  	};
};

// P:x,y->(<x|P_j y>)_l
void skyrmion_constrain_adjucent(const real* restrict a, const real* restrict b, real* restrict r) {
	#pragma omp parallel for collapse(4)	
  	forall(u,x,y,z) {
	  	int i=INDEX(u,x,y,z);
	  	r[i]=dot3(a+3*i,b+3*i);
  	};
};

void skyrmion_middle(const real* restrict a, const real* restrict b, real* restrict r) {
	#pragma omp parallel for collapse(4)	
  	forall(u,x,y,z) {
	  	int i=INDEX(u,x,y,z);
	  	middle3(a+3*i,b+3*i,r+3*i);
  	};
};

// calculate point in between of b and c by interpolation of curve given by
// points a,b,c and d in that order
void skyrmion_middle_fourth_order(const real* restrict a, const real* restrict b, const real* restrict c, const real* restrict d, real* restrict r) {
	#pragma omp parallel for collapse(4)	
  	forall(u,x,y,z) {
	  	int i=INDEX(u,x,y,z);
	  	middle_fourth_order3(a+3*i,b+3*i,c+3*i,d+3*i,r+3*i);
  	};
};

void skyrmion_middle_third_order(const real* restrict a, const real* restrict b, const real* restrict c, real* restrict r) {
	#pragma omp parallel for collapse(4)	
  	forall(u,x,y,z) {
	  	int i=INDEX(u,x,y,z);
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
		add_random_vector(noise*(m-n),size,p+k*size);
		normalize(p+k*size);
	};
	// Заполняем пробелы
	skyrmion_geodesic_rec(noise,p,n,k);
	skyrmion_geodesic_rec(noise,p,k,m);
}

void skyrmion_geodesic(real noise, int sizep, real* p) { 
	skyrmion_geodesic_rec(noise, p, 0, sizep-1); 
}

void two_point_tangent0(const real* restrict a, const real* restrict b, real* restrict r) {
	#pragma omp parallel for collapse(4)	
	forall(u,x,y,z) {	
		int i=INDEX(u,x,y,z)*3;
		sub3(b+i,a+i,r+i); 
		tangent3(a+i,r+i); 
		//normalize3(r+i);
	};
};

void two_point_tangent1(const real* restrict a, const real* restrict b, real* restrict r) {
	#pragma omp parallel for collapse(4)	
	forall(u,x,y,z) {	
		int i=INDEX(u,x,y,z)*3;
		sub3(b+i,a+i,r+i); 
		tangent3(b+i,r+i); 
		//normalize3(r+i);
	};
};

// tangent r to path defined by three consequative points a,b,c
void three_point_tangent(const real* restrict a, const real* restrict b, const real* restrict c, real* restrict r) {
	#pragma omp parallel for collapse(4)	
	forall(u,x,y,z) {	
		int i=INDEX(u,x,y,z)*3;
		sub3(c+i,a+i,r+i); 
		tangent3(b+i,r+i); 
		//normalize3(r+i);
	};
};

void three_point_tangent_stable(real ea, real eb, real ec, const real* restrict a, const real* restrict b, const real* restrict c, real* restrict r) {
	if(ea<eb && eb<ec) {
		#pragma omp parallel for collapse(4)	
		forall(u,x,y,z) {	
			int i=INDEX(u,x,y,z)*3;
			sub3(c+i,b+i,r+i); 
			tangent3(b+i,r+i); 
		};
	} else if(ea>eb && eb>ec) {
		#pragma omp parallel for collapse(4)	
		forall(u,x,y,z) {	
			int i=INDEX(u,x,y,z)*3;
			sub3(b+i,a+i,r+i); 
			tangent3(b+i,r+i); 
		};
	} else {
		real w1=rabs(eb-ea); real w2=rabs(ec-eb);
		if(w1>w2) { real t=w1; w1=w2; w2=t; }; // w1=min(|eb-ea|,|ec-eb|), w2=max(...)
		if(ec<ea) { real t=w1; w1=w2; w2=t; }; // w1=min if ec>ea  w1=max otherwise
		#pragma omp parallel for collapse(4)	
		forall(u,x,y,z) {	
			int i=INDEX(u,x,y,z)*3;
			real t1[3],t2[3];
			sub3(b+i,a+i,t1); 
			sub3(c+i,b+i,t2);
			for3(j) (r+i)[j]=w1*t1[j]+w2*t2[j];
			tangent3(b+i,r+i); 
		};
	}
};

// tangent r to path defined by three consequative points a,b,c
void three_point_tangent_mean(const real* restrict a, const real* restrict b, const real* restrict c, real* restrict r) {
	#pragma omp parallel for collapse(4)	
	forall(u,x,y,z) {	
		real t1[3],t2[3];
		int i=INDEX(u,x,y,z)*3;
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
void three_point_equalize(const real* restrict a, const real* restrict b, real* restrict r) {
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

void three_point_equalizer(const real* restrict a, const real* restrict c, const real* restrict b, real* restrict r) {
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
void append_skyrmion(const real center[3], real distance, real winding, 
	real rotation, real* restrict data) 
{
	real field[3]={0,0,1}; 
	#pragma omp parallel for collapse(4)	
	forall(u,x,y,z) {	
		int i=INDEX(u,x,y,z)*3;
		real vec[3]; COORDS(u,x,y,z,vec);
		sub3(vec,center,vec);
		//real elevation=dot3(vec,magnetic_field)/hnorm;
		//if(rabs(elevation)<distance) continue;
		//mult_minus3(elevation, magnetic_field, vec);
		real dist=rsqrt(normsq3(vec));
		if(dist>distance) continue;
		if(dist==0) {
			mult3(-1,data+i,data+i);
			continue;
		};
		multinv3(dist, vec, vec);
		// First rotation is around vec
		dist/=distance; dist=1-dist; dist*=M_PI_2;
		real sinalpha, cosalpha; 
		rsincos(dist*winding,&sinalpha,&cosalpha); 
		real q1[4]; q1[0]=cosalpha; mult3(sinalpha,vec,q1+1);
		// Second rotation is in the plane containing vec and field
		rsincos(dist*rotation,&sinalpha,&cosalpha); 
		real q2[4]; q2[0]=cosalpha; cross3(field,vec,q2+1);
		mult3(sinalpha,q2+1,q2+1);
		// combined rotation
		real q[4]; quaternion_product(q2,q1,q);
		// applying rotation
		q1[0]=0; copy3(data+i, q1+1);
		quaternion_product(q,q1,q2);
		q[1]=-q[1]; q[2]=-q[2]; q[3]=-q[3];
		quaternion_product(q2,q,q1);
		//fprintf(stderr,"Error %"RF"g\n",rsqrt(normsq3(q1+1))-1);
		//assert(rabs(rsqrt(normsq3(q1+1))-1)<1e-8);
		copy3(q1+1,data+i);
	};
}


// axis: 0=x, 1=y, 2=z
void group_generator(const real* restrict spins, 
int axis, real* restrict gen) {
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
		fprintf(stderr, COLOR_BOLD COLOR_RED"Error:"COLOR_RESET 
			"Axis %d has even length %d\n",axis,N);
		exit(1);
	};
// computing kernel	
	real* kernel=malloc(sizeof(real)*N); assert(kernel);
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



void fourier_table(const real* restrict angles, real* restrict table) {
	#pragma omp parallel for collapse(4)	
	forall(u,x,y,z) {
		int i=INDEX(u,x,y,z);
		rsincos(angles[2*i+0],table+4*i+0,table+4*i+1);
		rsincos(angles[2*i+1],table+4*i+2,table+4*i+3);
	};
};

// transform angles to vector on sphere
// (x,y,z)=R(phi,theta)
void angles_to_vector(const real* restrict table, const real* restrict angles, real* restrict vectors) {
	#pragma omp parallel for collapse(4)	
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