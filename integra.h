#ifndef INTEGRA_H
#define INTEGRA_H

#include "vector.h"

void runge_kutta(
	int N, 
	void (*F)(real* x, real* g),
	real T,
	real* X
	);

#endif