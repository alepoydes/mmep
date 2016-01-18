#ifndef OCTAVE_H
#define OCTAVE_H

#include "vector.h"
#include <stdio.h>

void oct_save_hessian(FILE* file);

void oct_save_linear(FILE* file);

void oct_save_state(FILE* file, char* name, real* data);

void oct_save_vertices(FILE* file);

void oct_save_init(FILE* file);

void oct_save_finish(FILE* file);

#endif