#ifndef INTEGRA_H
#define INTEGRA_H

#include "vector.h"

void runge_kutta(
	int N, 
	void (*F)(real* x, real* g),
	real T,
	real* X
	);

real runge_kutta_implicit(
	int N, 
	void (*F)(real* x, real* g),
	real T,
	int D,
	const real* restrict A,
	const real* restrict B,
	real* X,
	int tol,
	int max_iter
	);

real radau_integrator(
	int N, 
	void (*F)(real* x, real* g),
	real T,
	real* X,
	int tol,
	int max_iter
	);

real gauss_integrator(
	int N, 
	void (*F)(real* x, real* g),
	real T,
	real* X,
	int tol,
	int max_iter
	);

#endif