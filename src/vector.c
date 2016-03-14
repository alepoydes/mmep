#include "vector.h"

#include <stdlib.h>
#include <math.h>

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
