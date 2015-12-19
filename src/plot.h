#include <stdio.h>

#include "vector.h"

void plot_field3(FILE* file, const real* restrict a);
void plot_path(FILE* file, int sizep, const real* restrict mep);
void animate_path(FILE* file, int sizep, const real* restrict mep);
void plot_bounds(real bounds[3][2]);