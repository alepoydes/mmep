#include "vector.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

real normsq(int n, const real* a) {
	real nrm=0;
	#pragma omp parallel for reduction(+:nrm)
	for(int k=0; k<n; k++) nrm+=a[k]*a[k];
	return nrm;
};

real distsq(int n, const real* a, const real* b) {
	real nrm=0;
	#pragma omp parallel for reduction(+:nrm)
	for(int k=0; k<n; k++) nrm+=(a[k]-b[k])*(a[k]-b[k]);
	return nrm;
};

// n is three times smaller than for distsq
real dist_sphere_sq(int n, const real* a, const real* b) {
	real nrm=0;
	#pragma omp parallel for reduction(+:nrm)
	for(int k=0; k<n; k++) nrm+=dist_sphere_sq3(a+3*k,b+3*k);
	return nrm;	
};

real dot(int n, const real* a, const real* b) {
	real nrm=0;
	#pragma omp parallel for reduction(+:nrm)
	for(int k=0; k<n; k++) nrm+=a[k]*b[k];
	return nrm;	
};

void mult_sub(int n, real a, const real* restrict b, real* restrict c) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) c[k]-=a*b[k];
};
void mult_add(int n, real a, const real* restrict b, real* restrict c) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) c[k]+=a*b[k];
};

void random_vector(int n, real* a) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) a[k]=(real)rand()/(real)(RAND_MAX)-0.5;
};

void add_random_vector(real alpha, int n, real* a) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) a[k]+=alpha*((real)rand()/(real)(RAND_MAX)-0.5);	
};

void zero_vector(int n, real* a) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) a[k]=0;
};

void sub_const(int n, real a, real* c) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) c[k]-=a;
};

void copy_vector(int n, const real* restrict a, real* restrict b) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) b[k]=a[k];
};

void add_mult(int n, const real* restrict a, real b, real* restrict c) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) c[k]=b*c[k]+a[k];
}; 

void sub_inplace(int n, const real* restrict a, real* restrict c) {
	#pragma omp parallel for
	for(int k=0; k<n; k++) c[k]-=a[k];
}; 
void sub(int n, const real* restrict a, const real* restrict b, real* restrict c) {
	#pragma omp parallel for
	for(int k=0; k<n; k++) c[k]=a[k]-b[k];
};

void add_inplace(int n, const real* restrict a, real* restrict c) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) c[k]+=a[k];
}; 

void negate_inplace(int n, real* restrict c) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) c[k]=-c[k];
}; 

void negate_div(int n, real* restrict a, real* restrict c) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) 
		if(a[k]!=0) c[k]/=-a[k]; else c[k]=0;
}; 

void add_constant_inplace(int n, real a, real* restrict c) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) c[k]+=a;
};

void linear_comb(int n, real a, const real* restrict b, real c, real* restrict d, real* restrict e) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) e[k]=a*b[k]+c*d[k];
}; 

void const_div_inplace(int n, real a, real* restrict c) {
	#pragma omp parallel for 
	for(int k=0; k<n; k++) c[k]/=a;
};

// Orthogonalize family of K vectors V each of the length N.
// Vectors are not normalized.
// The resulting V[k] is the projection of input V[k]
// onto orthogonal subspace to V[j] for all j<k.
void gram_schmidt(int N, int K, real* restrict* V) {
	real* sqnorms=alloca(sizeof(real)*K);
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