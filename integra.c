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
	void (*F)(real* x, real* g),
	real T,
	real* X
	)
{
	// Allocate buffers
	real* k=malloc(sizeof(real)*N); assert(k);
	real* y=malloc(sizeof(real)*N); assert(y);
	real* g=malloc(sizeof(real)*N); assert(g);
	// Prepare data
	copy_vector(N,X,y);
	// first
	F(X,g);
	mult_add(N,T/6,g,y);
	// second
	copy_vector(N,X,k); 
	mult_add(N,T/2,g,k);
	F(k,g);
	mult_add(N,T/3,g,y);
	// third
	copy_vector(N,X,k); 
	mult_add(N,T/2,g,k);
	F(k,g);
	mult_add(N,T/3,g,y);
	// forth
	copy_vector(N,X,k); 
	mult_add(N,T,g,k);
	F(k,g);
	mult_add(N,T/6,g,y);
	// Copy result
	copy_vector(N,y,X);
	// Deallocate buffers
	free(k); free(y); free(g);
};