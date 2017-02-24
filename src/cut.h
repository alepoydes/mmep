#ifndef CUT_H
#define CUT_H

#include "vector.h"

#include <stdio.h>

typedef enum {
    cut_plane, cut_sphere, cut_complement, cut_union, cut_intersection, cut_universe
} cut_t;

typedef struct primitive { 
    cut_t type;
    union {
        struct {
            real normal[3];
            real shift;
        } plane;
        struct {
            real center[3];
            real radius;
        } sphere;
        struct {
            struct primitive* a;
        } complement;
        struct {
            struct primitive* a;
            struct primitive* b;
        } setunion;
        struct {
            struct primitive* a;
            struct primitive* b;
        } intersection;
    } data;
} primitive_t;

primitive_t* new_plane(real point[3], real normal[3]);
primitive_t* new_sphere(real center[3], real radius);
primitive_t* new_complement(primitive_t* a);
primitive_t* new_union(primitive_t* a, primitive_t* b);
primitive_t* new_intersection(primitive_t* a, primitive_t* b);
primitive_t* new_universe();

bool does_belong_to_primitive(primitive_t* primitive, real point[3]);
void fprint_primitive(FILE* file, primitive_t* pr);

#endif