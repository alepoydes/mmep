#ifndef VECTOR_H
#define VECTOR_H

#include <math.h>
#include <float.h>

#ifdef DOUBLE
	#define real double
	#define rsqrt(x) sqrt(x)
	#define rabs(x) fabs(x)
	#define rsincos(x,y,z) sincos(x,y,z)
	#define RF "l"
	#define EPSILON DBL_EPSILON
#else
#ifdef QUAD
	#define real long double
	#define rsqrt(x) sqrtl(x)
	#define rabs(x) fabsl(x)
	#define rsincos(x,y,z) sincosl(x,y,z)
	#define RF "L"
	#define EPSILON LDBL_EPSILON
#else 
	#define real float
	#define rsqrt(x) sqrtf(x)
	#define rabs(x) fabsf(x)
	#define rsincos(x,y,z) sincosf(x,y,z)
	#define RF ""
	#define EPSILON FLT_EPSILON
#endif
#endif

// считает длину вектора
#define normsq3(x) ((x)[0]*(x)[0]+(x)[1]*(x)[1]+(x)[2]*(x)[2])
// считает скалярное произведение векторов
#define dot3(x,y) (x)[0]*(y)[0]+(x)[1]*(y)[1]+(x)[2]*(y)[2]
//#define normalize3(x) { real t=(3-normsq3(x))*0.5; (x)[0]*=t; (x)[1]*=t; (x)[2]*=t; }
//#define normalize3(x) { real t=(1+normsq3(x))/2; (x)[0]/=t; (x)[1]/=t; (x)[2]/=t; }
#define normalize3(x) { real t=rsqrt(normsq3(x)); if(t>0) { (x)[0]/=t; (x)[1]/=t; (x)[2]/=t; }; }
#define middle3(x,y,z) { (z)[0]=((x)[0]+(y)[0])/2; (z)[1]=((x)[1]+(y)[1])/2; (z)[2]=((x)[2]+(y)[2])/2; normalize3(z); }
#define sub3(x,y,z) { (z)[0]=(x)[0]-(y)[0]; (z)[1]=(x)[1]-(y)[1]; (z)[2]=(x)[2]-(y)[2]; }
// проецирует y на касательное подпространство к x
#define tangent3(x,y) {	real t=dot3(x,y); (y)[0]-=t*(x)[0]; (y)[1]-=t*(x)[1]; (y)[2]-=t*(x)[2]; }
#define cross_minus3(a,b,c) { (c)[0]-=(a)[1]*(b)[2]-(a)[2]*(b)[1]; (c)[1]-=(a)[2]*(b)[0]-(a)[0]*(b)[2]; (c)[2]-=(a)[0]*(b)[1]-(a)[1]*(b)[0]; }
#define cross_plus3(a,b,c) { (c)[0]+=(a)[1]*(b)[2]-(a)[2]*(b)[1]; (c)[1]+=(a)[2]*(b)[0]-(a)[0]*(b)[2]; (c)[2]+=(a)[0]*(b)[1]-(a)[1]*(b)[0]; }
#define mult_minus3(a,b,c) { (c)[0]-=(a)*(b)[0]; (c)[1]-=(a)*(b)[1]; (c)[2]-=(a)*(b)[2]; }
#define mult_plus3(a,b,c) { (c)[0]+=(a)*(b)[0]; (c)[1]+=(a)*(b)[1]; (c)[2]+=(a)*(b)[2]; }
#define copy3(a,c) { (c)[0]=(a)[0]; (c)[1]=(a)[1]; (c)[2]=(a)[2]; }

real distsq(int n, const real* a, const real* b);
real normsq(int n, const real* a);
real dot(int n, const real* a, const real* b);
void mult_sub(int n, real a, const real* b, real* c);
void mult_add(int n, real a, const real* b, real* c);
void add_mult(int n, const real* a, real b, real* c); 
void sub_const(int n, real a, real* c);
void random_vector(int n, real* a);
void zero_vector(int n, real* a);
void copy_vector(int n, const real* a, real* b);
void sub_inplace(int n, const real* restrict a, real* restrict c);
void add_inplace(int n, const real* restrict a, real* restrict c);
void negate_inplace(int n, real* restrict c);

#endif