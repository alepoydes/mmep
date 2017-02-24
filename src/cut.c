#include "cut.h"

#include <stdlib.h>
#include <stdio.h>

primitive_t* new_plane(real point[3], real normal[3]) {
    primitive_t* result=(primitive_t*)malloc(sizeof(primitive_t));
    result->type=cut_plane;
    let3(result->data.plane.normal, normal);
    result->data.plane.shift=-normal[0]*point[0]-normal[1]*point[1]-normal[2]*point[2];
    return result;
};

primitive_t* new_sphere(real center[3], real radius) {
    primitive_t* result=(primitive_t*)malloc(sizeof(primitive_t));
    result->type=cut_sphere;
    let3(result->data.sphere.center, center);
    result->data.sphere.radius=radius;
    return result;
};

primitive_t* new_complement(primitive_t* a) {
    primitive_t* result=(primitive_t*)malloc(sizeof(primitive_t));
    result->type=cut_complement;
    result->data.complement.a=a;
    return result;
};

primitive_t* new_union(primitive_t* a, primitive_t* b) {
    primitive_t* result=(primitive_t*)malloc(sizeof(primitive_t));
    result->type=cut_union;
    result->data.setunion.a=a;
    result->data.setunion.b=b;
    return result;
};

primitive_t* new_intersection(primitive_t* a, primitive_t* b) {
    primitive_t* result=(primitive_t*)malloc(sizeof(primitive_t));
    result->type=cut_intersection;
    result->data.intersection.a=a;
    result->data.intersection.b=b;
    return result;
};

primitive_t* new_universe() {
    primitive_t* result=(primitive_t*)malloc(sizeof(primitive_t));
    result->type=cut_universe;
    return result;
};

bool does_belong_to_primitive(primitive_t* primitive, real point[3]) {
    real tmp[3];
    switch (primitive->type) {
        case cut_universe:
            return true;
        case cut_complement: 
            return !does_belong_to_primitive(primitive->data.complement.a, point);
        case cut_union: 
            return does_belong_to_primitive(primitive->data.setunion.a, point) 
                   || does_belong_to_primitive(primitive->data.setunion.b, point);
        case cut_intersection: 
            return does_belong_to_primitive(primitive->data.intersection.a, point) 
                   && does_belong_to_primitive(primitive->data.intersection.b, point);
        case cut_sphere: 
            sub3(point, primitive->data.sphere.center, tmp);
            return rsqrt(normsq3(tmp))<=primitive->data.sphere.radius;
        case cut_plane: 
            return dot3(primitive->data.plane.normal, point)+primitive->data.plane.shift>=0;
        default:
            fprintf(stderr, "Error: unimplemented primitive\n");
            exit(1);
    };
};

void fprint_primitive(FILE* file, primitive_t* primitive) {
    switch (primitive->type) {
        case cut_universe:
            fprintf(file, "universe");
            break;
        case cut_complement: 
            fprintf(file, "!");
            fprint_primitive(file,primitive->data.complement.a);
            break;
        case cut_union: 
            fprintf(file, "(");
            fprint_primitive(file,primitive->data.setunion.a);
            fprintf(file, ")+(");
            fprint_primitive(file,primitive->data.setunion.b);
            fprintf(file, ")");
            break;
        case cut_intersection: 
            fprintf(file, "(");
            fprint_primitive(file,primitive->data.intersection.a);
            fprintf(file, ")*(");
            fprint_primitive(file,primitive->data.intersection.b);
            fprintf(file, ")");
            break;
        case cut_sphere: 
            fprintf(file, "sphere {%"RF"g, %"RF"g, %"RF"g} %"RF"g", 
                RT(primitive->data.sphere.center[0]),RT(primitive->data.sphere.center[1]),RT(primitive->data.sphere.center[2]),
                RT(primitive->data.sphere.radius));
            break;
        case cut_plane: 
            fprintf(file, "plane %"RF"g {%"RF"g, %"RF"g, %"RF"g}", 
                RT(primitive->data.plane.shift),
                RT(primitive->data.plane.normal[0]),RT(primitive->data.plane.normal[1]),RT(primitive->data.plane.normal[2]));
            break;
        default:
            fprintf(stderr, "Error: unimplemented primitive\n");
            exit(1);
    };
};
