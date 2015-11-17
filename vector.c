#include <stdlib.h>
#include <math.h>

#include "vector.h"

real normsq(int n, const real* a) {
	real nrm=0;
	for(int k=0; k<n; k++) nrm+=a[k]*a[k];
	return nrm;
};

real dot(int n, const real* a, const real* b) {
	real nrm=0;
	for(int k=0; k<n; k++) nrm+=a[k]*b[k];
	return nrm;	
};

void mult_sub(int n, real a, const real* restrict b, real* restrict c) {
	for(int k=0; k<n; k++) c[k]-=a*b[k];
};
void mult_add(int n, real a, const real* restrict b, real* restrict c) {
	for(int k=0; k<n; k++) c[k]+=a*b[k];
};

void random_vector(int n, real* a) {
	for(int k=0; k<n; k++) a[k]=(real)rand()/(real)(RAND_MAX)-0.5;
};

void zero_vector(int n, real* a) {
	for(int k=0; k<n; k++) a[k]=0;
};

void sub_const(int n, real a, real* c) {
	for(int k=0; k<n; k++) c[k]-=a;
};

void copy_vector(int n, const real* restrict a, real* restrict b) {
	for(int k=0; k<n; k++) b[k]=a[k];
};

void add_mult(int n, const real* restrict a, real b, real* restrict c) {
	for(int k=0; k<n; k++) c[k]=b*c[k]+a[k];
}; 

void sub_inplace(int n, const real* restrict a, real* restrict c) {
	for(int k=0; k<n; k++) c[k]-=a[k];
}; 

void negate_inplace(int n, real* restrict c) {
	for(int k=0; k<n; k++) c[k]=-c[k];
}; 

void vector_copy(int n, const real* a, real* b) {
	for(int k=0; k<n; k++) b[k]=a[k];
}