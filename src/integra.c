#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "vector.h"
#include "debug.h"
#include "integra.h"

/*
Solve ODE doing Runge-Kutta
ODE is defined by F(X):
	X'(t)=F(X(t))
Arguments:
	N - rank of X
	F(X,G) - compute G=X' given X,
	T - time step,
	X=X(0) - initial state,
Result X(T) is stored in X.
*/
void runge_kutta(
	int N, 
	void (*F)(const real* x, real* g, realp* E),
	real T,
	real* X, 
	realp* E,
	int* iter
	)
{
	// Allocate buffers
	real* k=malloc(sizeof(real)*N); assert(k);
	real* y=malloc(sizeof(real)*N); assert(y);
	real* g=malloc(sizeof(real)*N); assert(g);
	// Prepare data
	copy_vector(N,X,y);
	// first
	F(X,g,E);
	mult_add(N,T/6,g,y);
	// second
	copy_vector(N,X,k); 
	mult_add(N,T/2,g,k);
	F(k,g,NULL);
	mult_add(N,T/3,g,y);
	// third
	copy_vector(N,X,k); 
	mult_add(N,T/2,g,k);
	F(k,g,NULL);
	mult_add(N,T/3,g,y);
	// forth
	copy_vector(N,X,k); 
	mult_add(N,T,g,k);
	F(k,g,NULL);
	mult_add(N,T/6,g,y);
	// Copy result
	copy_vector(N,y,X);
	if(iter) *iter=4;
	// Deallocate buffers
	free(k); free(y); free(g);
};

void euler(
	int N, 
	void (*F)(const real* x, real* g, realp* E),
	real T,
	real* X, 
	realp* E,
	int* iter
	)
{
	// Allocate buffers
	real* g=malloc(sizeof(real)*N); assert(g);
	// first
	F(X,g,E);
	mult_add(N,T,g,X);
	if(iter) *iter=1;
	// Deallocate buffers
	free(g); 
};


/* 
Do single step of implicit Runge-Kutta
for ODE 
	X'(t)=F(X(t))
with initial condition X(0),
the result is X(T).
Arguments:
	X - storage for X(0) and X(T)
	N - dim of X
	F(X,G) - compute vectorfield F(X)and store result in G,
	T - integration step,
	D - number of intermediate steps,
	A,B - Butcher tableau,
	A - matrix D x D
	B - vector of dim D
	tol - error tolerance
	max_iter - maximum number of iterations
*/

real runge_kutta_implicit(
	int N, 
	void (*F)(const real* x, real* g, realp* E),
	real T,
	int D,
	const real* restrict A,
	const real* restrict B,
	real* X,
	real tol,
	int max_iter,
	realp* E,
	int* iter
	)
{
	// check arguments
	assert(D>0); assert(N>=0);
	// Allocate buffers
	real* mu=malloc(sizeof(real)*N*D); assert(mu);
	real* mu1=malloc(sizeof(real)*N*D); assert(mu1);
	real* g=malloc(sizeof(real)*N); assert(g);
	// Initial approximation
	real err=NAN;
	F(X, mu, E);
	if(iter) *iter=1;
	for(int d=1; d<D; d++) copy_vector(N,mu,mu+d*N);
	// Repeat until converge
	for(int i=0; i<max_iter; i++) {
		// compute next iteration
		for(int d=0; d<D; d++) {
			copy_vector(N,X,g);
			for(int j=0; j<D; j++)
				mult_add(N,T*A[d*D+j],mu+j*N,g);
			F(g,mu1+d*N,NULL);
			if(iter) (*iter)++;
		};
		// swap buffers
		real* tmp=mu1; mu1=mu; mu=tmp;
		// check if converged
		err=rpsqrt(distsq(N*D, mu, mu1));
		//fprintf(stderr, "RK: iter %d err %"RF"g\n",i,err);
		if(err<tol) break;
	}; // intermediate steps are in mu
	free(mu1); free(g);
	// compute final step
	for(int d=0; d<D; d++) mult_add(N,T*B[d],mu+d*N,X);
	// deallocating
	free(mu); 
	return err;
};

real radau_integrator(
	int N, 
	void (*F)(const real* x, real* g, realp* E),
	real T,
	real* X,
	int tol,
	int max_iter, 
	realp* E,
	int* iter
	)
{
	real r=rsqrt(6);
	real B[3]={(16-r)/36,(16+r)/36,1./9};
	real A[9]={(16.-r)/72, (328-167*r)/1800,(-2+3*r)/450,
		(328+167*r)/1800,(16+r)/72,(-2-3*r)/450,
		(85-10*r)/180,(85+10*r)/180,1./18};
	return runge_kutta_implicit(N, F, T, 3, A, B, X, tol, max_iter, E, iter);
};

real gauss_integrator(
	int N, 
	void (*F)(const real* x, real* g, realp* E),
	real T,
	real* X,
	int tol,
	int max_iter, 
	realp* E,
	int* iter
	)
{
	real r=rsqrt(3);
	real B[2]={0.5,0.5};
	real A[4]={1./4, 1./4-r/6, 1./4+r/6, 1./4};
	return runge_kutta_implicit(N, F, T, 2, A, B, X, tol, max_iter, E, iter);
};