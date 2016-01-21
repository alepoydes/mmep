#define _GNU_SOURCE
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "vector.h"
#include "skyrmion.h"

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
// Dzyaloshinskii Moriya vector for every pair of atoms
real* dzyaloshinskii_moriya_vector=NULL;
real* initial_state=NULL;
real* final_state=NULL;

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
	#pragma omp parallel for collapse(4)
	forall(u,x,y,z) for3(j) out[INDEX(u,x,y,z)*3+j]=magnetic_field[j];
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
	
void skyrmion_geodesic_rec(real* p, int n, int m) {
	int size=SIZE*3;
	if(n+1>=m) return;
	int k=(n+m)/2; // Точка разбиения
	// Находим середину на прямой между n и m и проецируем на сферы
	skyrmion_middle(p+n*size,p+m*size,p+k*size);
	// Заполняем пробелы
	skyrmion_geodesic_rec(p,n,k);
	skyrmion_geodesic_rec(p,k,m);
}

void skyrmion_geodesic(int sizep, real* p) { 
	skyrmion_geodesic_rec(p, 0, sizep-1); 
}

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
void append_skyrmion(const real center[3], real distance, int winding, 
	int rotation, real* restrict data) 
{
	real field[3]; copy3(magnetic_field, field); normalize3(field);
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
			if((rotation+winding)%2!=0) mult3(-1,data+i,data+i);
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
		real q[4]; quaternion_product(q1,q2,q);
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