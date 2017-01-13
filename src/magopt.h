#ifndef MAGOPT_H
#define MAGOPT_H

#include "vector.h"

#include <stdio.h>

//////////////////////////////////////////////////////
// Single image optimization

int skyrmion_steepest_descent(real* __restrict__ x, 
int mode, real mode_param, real epsilon, int max_iter);

int skyrmion_better_steepest_descent(
real* __restrict__ x, real epsilon, int max_iter
);

int skyrmion_minimize(real* __restrict__ x, real epsilon, int max_iter);
//////////////////////////////////////////////////////
// MEP optimization

extern int sizep; // Number of nodes on path
extern int debug_plot_path;
extern int remove_zero_modes;
extern int flat_distance;
extern int post_optimization;
extern int single_maximum;

extern realp *distance;
extern realp *energy;
extern realp *diff;
extern realp *tdiff;
extern realp *inflation;

void path_steepest_descent_init(int max_sizep);
void path_steepest_descent_deinit();

int path_steepest_descent(real* __restrict__ path, int mode, 
	real mode_param, real epsilon, int max_iter);
void energy_evaluate(real* path);
void energy_display(FILE* file);

#endif // MAGOPT_H