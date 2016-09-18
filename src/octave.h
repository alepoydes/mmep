#ifndef OCTAVE_H
#define OCTAVE_H

#include "vector.h"
#include <stdio.h>

void oct_save_hessian(FILE* file);

void oct_save_linear(FILE* file);

void oct_save_state(FILE* file, const char* name, real* data);

void oct_save_vertices(FILE* file);

void oct_save_init(FILE* file);

void oct_save_finish(FILE* file);

void oct_save_real(FILE* file, const char* name, real value);
void oct_save_realp(FILE* file, const char* name, realp value);

void oct_save_field(FILE* file, const char* name, real* data);

void oct_save_path(FILE* file, const char* name, real* data, int sizep);

void oct_save_vector(FILE* file, const char* name, real* data, int length);
void oct_save_vectorp(FILE* file, const char* name, realp* data, int length);
void oct_save_vector_int(FILE* file, const char* name, int* data, int length);

void oct_save_matrix(FILE* file, const char* name, real* data, int height, int width);
void oct_save_matrixp(FILE* file, const char* name, realp* data, int height, int width);
void oct_save_matrix_int(FILE* file, const char* name, int* data, int height, int width);

#endif