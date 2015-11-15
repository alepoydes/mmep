#include <stdlib.h>
#include <math.h>

#include "vector.h"

float normsq(int n, const float* a) {
	float nrm=0;
	for(int k=0; k<n; k++) nrm+=a[k]*a[k];
	return nrm;
};

float dot(int n, const float* a, const float* b) {
	float nrm=0;
	for(int k=0; k<n; k++) nrm+=a[k]*b[k];
	return nrm;	
};

void mult_sub(int n, float a, const float* restrict b, float* restrict c) {
	for(int k=0; k<n; k++) c[k]-=a*b[k];
};
void mult_add(int n, float a, const float* restrict b, float* restrict c) {
	for(int k=0; k<n; k++) c[k]+=a*b[k];
};

void random_vector(int n, float* a) {
	for(int k=0; k<n; k++) a[k]=(float)rand()/(double)(RAND_MAX)-0.5;
};

void zero_vector(int n, float* a) {
	for(int k=0; k<n; k++) a[k]=0;
};

void sub_const(int n, float a, float* c) {
	for(int k=0; k<n; k++) c[k]-=a;
};

void copy_vector(int n, const float* restrict a, float* restrict b) {
	for(int k=0; k<n; k++) b[k]=a[k];
};

void add_mult(int n, const float* restrict a, float b, float* restrict c) {
	for(int k=0; k<n; k++) c[k]=b*c[k]+a[k];
}; 

void sub_inplace(int n, const float* restrict a, float* restrict c) {
	for(int k=0; k<n; k++) c[k]-=a[k];
}; 

void negate_inplace(int n, float* restrict c) {
	for(int k=0; k<n; k++) c[k]=-c[k];
}; 