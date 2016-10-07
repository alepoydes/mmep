#include "vector.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

real NORMEPS2=0;//1e-14;

real* ralloc(int n) {
	real* buf=(real*)malloc(sizeof(real)*n); 
	assert(buf);
	return buf;
};

real normalize3(real* x) { 
	real t=normsq3(x); 
	if(t==0) return 0;
	real e=rabs(t-1);
	if (e>NORMEPS2) { 
		real f=rsqrt(t); 
		(x)[0]/=f; (x)[1]/=f; (x)[2]/=f; 
		return 0;
	}; 
	return e;//EPSILON;//rabs(normsq3(x)-1);
};

real dot3(const real* __restrict__ x, const real* __restrict__ y) {
	return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
}; 

void cross3(const real* __restrict__ a,const real* __restrict__ b,real* __restrict__ c) { 
	(c)[0]=(a)[1]*(b)[2]-(a)[2]*(b)[1]; 
	(c)[1]=(a)[2]*(b)[0]-(a)[0]*(b)[2]; 
	(c)[2]=(a)[0]*(b)[1]-(a)[1]*(b)[0]; 
};

void invertmatrix3(const real* mat, real* inv) {
	real comp[3][3];
	for3(j) for3(k) 
		comp[j][k]=mat[3*((j+1)%3)+(k+1)%3]*mat[3*((j+2)%3)+(k+2)%3]-mat[3*((j+1)%3)+(k+2)%3]*mat[3*((j+2)%3)+(k+1)%3];
	real det=0; for3(j) det+=mat[j]*comp[0][j];
	for3(j) for3(k) inv[3*j+k]=comp[k][j]/det;
};

void matrixmult3(const real* a, const real* b, real* prod) {
	for3(j) for3(k) {
		real sum=0; for3(n) sum+=a[3*j+n]*b[3*n+k];
		prod[3*j+k]=sum;
	};
};

realp normsq(int n, const real* a) {
	realp nrm=0;
	#pragma omp parallel for reduction(+:nrm)
	for(int k=0; k<n; k++) nrm+=a[k]*a[k];
	return nrm;
};

realp distsq(int n, const real* a, const real* b) {
	realp nrm=0;
	#pragma omp parallel for reduction(+:nrm)
	for(int k=0; k<n; k++) nrm+=(a[k]-b[k])*(a[k]-b[k]);
	return nrm;
};

// n is three times smaller than for distsq
realp dist_sphere_sq(int n, const real* a, const real* b) {
	realp nrm=0;
	#pragma omp parallel for reduction(+:nrm)
	for(int k=0; k<n; k++) nrm+=dist_sphere_sq3(a+3*k,b+3*k);
	return nrm;	
};

real dist_sphere_sq3(const real* a,const real* b) {
	real p=dot3(a,b);
	if(p>1) p=1; else if(p<-1) p=-1;
	real d=racos(p);
	return d*d;
}; 

realp dot(int n, const real* a, const real* b) {
	realp nrm=0;
	#pragma omp parallel for reduction(+:nrm)
	for(int k=0; k<n; k++) nrm+=a[k]*b[k];
	return nrm;	
};

void mult_sub(int n, real a, const real* __restrict__ b, real* __restrict__ c) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) c[k]-=a*b[k];
};

void mult_sub_ext(int n, real a, const real* b, const real* c, real* d) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) d[k]=c[k]-a*b[k];
};

void mult_add(int n, real a, const real* __restrict__ b, real* __restrict__ c) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) c[k]+=a*b[k];
};

void random_vector(int n, real* a) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) a[k]=random_real()-0.5;
};

void add_random_vector(real alpha, int n, const real* a, real* b) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) b[k]=a[k]+alpha*(random_real()-0.5);	
};

// a and b are matrices [n x 3]
void add_random_cone(real alpha, int n, const real* a, real* b) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) {
		const real* A=a+3*k; real* B=b+3*k;
		real phi=random_real()*2*M_PI;
		real z=random_real()*2-1;
		real zc=rsqrt(1-z*z);
		real cphi, sphi; rsincos(phi, &cphi, &sphi);
		B[0]=A[0]+alpha*cphi*zc;
		B[1]=A[1]+alpha*sphi*zc;
		B[2]=A[2]+alpha*z;
	};
};

void add_random_cone3(real alpha, const real* A, real* B) {
	real phi=random_real()*2*M_PI;
	real z=random_real()*2-1;
	real zc=rsqrt(1-z*z);
	real cphi, sphi; rsincos(phi, &cphi, &sphi);
	B[0]=A[0]+alpha*cphi*zc;
	B[1]=A[1]+alpha*sphi*zc;
	B[2]=A[2]+alpha*z;
};

void zero_vector(int n, real* a) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) a[k]=0;
};

void sub_const(int n, real a, real* c) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) c[k]-=a;
};

void copy_vector(int n, const real* __restrict__ a, real* __restrict__ b) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) b[k]=a[k];
};

void add_mult(int n, const real* __restrict__ a, real b, real* __restrict__ c) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) c[k]=b*c[k]+a[k];
}; 

void sub_inplace(int n, const real* __restrict__ a, real* __restrict__ c) {
	#pragma omp parallel for
	for(int k=0; k<n; k++) c[k]-=a[k];
}; 
void sub(int n, const real* __restrict__ a, const real* __restrict__ b, real* __restrict__ c) {
	#pragma omp parallel for
	for(int k=0; k<n; k++) c[k]=a[k]-b[k];
};

void add_inplace(int n, const real* __restrict__ a, real* __restrict__ c) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) c[k]+=a[k];
}; 

void negate_inplace(int n, real* __restrict__ c) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) c[k]=-c[k];
}; 

void negate_div(int n, real* __restrict__ a, real* __restrict__ c) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) 
		if(a[k]!=0) c[k]/=-a[k]; else c[k]=0;
}; 

void add_constant_inplace(int n, real a, real* __restrict__ c) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) c[k]+=a;
};

void linear_comb(int n, real a, const real* __restrict__ b, real c, real* __restrict__ d, real* __restrict__ e) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) e[k]=a*b[k]+c*d[k];
}; 

void const_div_inplace(int n, real a, real* __restrict__ c) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) c[k]/=a;
};

// Orthogonalize family of K vectors V each of the length N.
// Vectors are not normalized.
// The resulting V[k] is the projection of input V[k]
// onto orthogonal subspace to V[j] for all j<k.
void gram_schmidt(int N, int K, real* __restrict__* V) {
	real sqnorms[K];
	for(int k=0; k<K; k++) {
		for(int j=0; j<k; j++) 
			if(sqnorms[j]>1e-7) {
			// subtract from V[k] projection onto V[j]
			real proj=dot(N,V[j],V[k]);
			mult_sub(N,proj/sqnorms[j],V[j],V[k]);
		};
		// Normalize V[k]
		sqnorms[k]=normsq(N,V[k]);
	};
}