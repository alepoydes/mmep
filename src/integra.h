#ifndef INTEGRA_H
#define INTEGRA_H

#include "vector.h"

void euler(
	int N, 
	void (*F)(const real* x, real* g, real* E),
	real T,
	real* X, 
	real* E,
	int* iter
	);

void runge_kutta(
	int N, 
	void (*F)(const real* x, real* g, real* E),
	real T,
	real* X, 
	real* E,
	int* iter
	);

real runge_kutta_implicit(
	int N, 
	void (*F)(const real* x, real* g, real* E),
	real T,
	int D,
	const real* restrict A,
	const real* restrict B,
	real* X,
	int tol,
	int max_iter, 
	real* E,
	int* iter
	);

real radau_integrator(
	int N, 
	void (*F)(const real* x, real* g, real* E),
	real T,
	real* X,
	int tol,
	int max_iter, 
	real* E,
	int* iter
	);

real gauss_integrator(
	int N, 
	void (*F)(const real* x, real* g, real* E),
	real T,
	real* X,
	int tol,
	int max_iter, 
	real* E,
	int* iter
	);

#endif