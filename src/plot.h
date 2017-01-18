#include <stdio.h>

#include "vector.h"

void plot_field3(FILE* file, const real* __restrict__ a);
void plot_path(FILE* file, int sizep, const real* __restrict__ mep);
void animate_path(FILE* file, int sizep, const real* __restrict__ mep);
void plot_bounds(real bounds[3][2]);
int load_path_from_gnuplot(FILE* file, real*(*allocate_image)());