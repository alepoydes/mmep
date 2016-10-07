#ifndef VECTOR_H
#define VECTOR_H

#define _GNU_SOURCE
#include <math.h>
#include <float.h>

#ifdef DOUBLE

	#define real double
	#define rsqrt(x) sqrt(x)
	#define rabs(x) fabs(x)
	#define racos(x) acos(x)
	#define rexp(x) exp(x)
#if defined(__clang__)
	#define rsincos(x,y,z) { *(y)=sin(x); *(z)=cos(x); }
#elif defined(__GNUC__) || defined(__GNUG__)
	#define rsincos(x,y,z) sincos(x,y,z)
#elif defined(_MSC_VER)
	#define rsincos(x,y,z) { *(y)=sin(x); *(z)=cos(x); }
#endif
	#define RF "l"
	#define SRF RF
	#define RT(x) (x)
	#define EPSILON DBL_EPSILON
	#define DIGITS DBL_DIG
	#define R_PI M_PI
	#define SPRINTF(...) sprintf(__VA_ARGS__)

#elif defined QUAD
	#include <quadmath.h>

	#define real __float128
	#define rsqrt(x) sqrtq(x)
	#define rabs(x) fabsq(x)
	#define racos(x) acosq(x)
	#define rexp(x) expq(x)
#if defined(__clang__)
	#define rsincos(x,y,z) { *(y)=sinq(x); *(z)=cosq(x); }
#elif defined(__GNUC__) || defined(__GNUG__)
	#define rsincos(x,y,z) sincosq(x,y,z)
#elif defined(_MSC_VER)
	#define rsincos(x,y,z) { *(y)=sinq(x); *(z)=cosq(x); }
#endif
	#define RF "L"
	#define SRF "Q"
	#define RT(x) (long double)(x)
	#define EPSILON FLT128_EPSILON
	#define DIGITS FLT128_DIG
	#define R_PI M_PIq
	#define SPRINTF(buf, ...) quadmath_snprintf(buf, sizeof(buf), __VA_ARGS__)

#else 

	#define real float
	#define rsqrt(x) sqrtf(x)
	#define rabs(x) fabsf(x)
	#define racos(x) acosf(x)
	#define rexp(x) expf(x)
#if defined(__clang__)
	#define rsincos(x,y,z) { *(y)=sinf(x); *(z)=cosf(x); }
#elif defined(__GNUC__) || defined(__GNUG__)
	#define rsincos(x,y,z) sincosf(x,y,z)
#elif defined(_MSC_VER)
	#define rsincos(x,y,z) { *(y)=sinf(x); *(z)=cosf(x); }
#endif
	#define RF ""
	#define SRF RF
	#define RT(x) (x)
	#define EPSILON FLT_EPSILON
	#define DIGITS FLT_DIG
	#define R_PI M_PI
	#define SPRINTF(...) sprintf(__VA_ARGS__)

#endif  

#define realp real
#define rpsqrt(x) rsqrt(x)
#define rpabs(x) rabs(x)
#define rpacos(x) acosl(x)
#define rpexp(x) expl(x)
#define rpsincos(x,y,z) rsincos(x,y,z)
#define RPF RF
#define PEPSILON EPSILON
#define PDIGITS DIGITS

real* ralloc(int n);
// считает длину вектора
#define normsq3(x) ((x)[0]*(x)[0]+(x)[1]*(x)[1]+(x)[2]*(x)[2])
// считает скалярное произведение векторов
real dot3(const real* x, const real* y);
extern real NORMEPS2;
real normalize3(real* x);
#define seminormalize3(factor,x) ({ real t=rsqrt(normsq3(x)); if(t>0) { real f=1-factor+factor/t; (x)[0]*=f; (x)[1]*=f; (x)[2]*=f; }; t>0?rabs((1-t)*(1-factor)):1; })
#define middle3(x,y,z) { (z)[0]=((x)[0]+(y)[0])/2; (z)[1]=((x)[1]+(y)[1])/2; (z)[2]=((x)[2]+(y)[2])/2; normalize3(z); }
#define middle_fourth_order3(a,b,c,d,z) { (z)[0]=(-(a)[0]+9*(b)[0]+9*(c)[0]-(d)[0])/16; (z)[1]=(-(a)[1]+9*(b)[1]+9*(c)[1]-(d)[1])/16; (z)[2]=(-(a)[2]+9*(b)[2]+9*(c)[2]-(d)[2])/16;  normalize3(z); }
#define middle_third_order3(a,b,c,z) { (z)[0]=(3*(a)[0]+6*(b)[0]-(c)[0])/8; (z)[1]=(3*(a)[1]+6*(b)[1]-(c)[1])/8; (z)[2]=(3*(a)[2]+6*(b)[2]-(c)[2])/8;  normalize3(z); }

#define sub3(x,y,z) { (z)[0]=(x)[0]-(y)[0]; (z)[1]=(x)[1]-(y)[1]; (z)[2]=(x)[2]-(y)[2]; }
// проецирует y на касательное подпространство к x
#define tangent3(x,y) {	real t=dot3(x,y); (y)[0]-=t*(x)[0]; (y)[1]-=t*(x)[1]; (y)[2]-=t*(x)[2]; }
#define reverse3(x,y) {	real t=1.5*dot3(x,y); (y)[0]-=t*(x)[0]; (y)[1]-=t*(x)[1]; (y)[2]-=t*(x)[2]; }
//#define cross3(a,b,c) { (c)[0]=(a)[1]*(b)[2]-(a)[2]*(b)[1]; (c)[1]=(a)[2]*(b)[0]-(a)[0]*(b)[2]; (c)[2]=(a)[0]*(b)[1]-(a)[1]*(b)[0]; }
void cross3(const real* a,const real* b,real* c);
#define cross_minus3(a,b,c) { (c)[0]-=(a)[1]*(b)[2]-(a)[2]*(b)[1]; (c)[1]-=(a)[2]*(b)[0]-(a)[0]*(b)[2]; (c)[2]-=(a)[0]*(b)[1]-(a)[1]*(b)[0]; }
#define cross_plus3(a,b,c) { (c)[0]+=(a)[1]*(b)[2]-(a)[2]*(b)[1]; (c)[1]+=(a)[2]*(b)[0]-(a)[0]*(b)[2]; (c)[2]+=(a)[0]*(b)[1]-(a)[1]*(b)[0]; }
#define mult_minus3(a,b,c) { (c)[0]-=(a)*(b)[0]; (c)[1]-=(a)*(b)[1]; (c)[2]-=(a)*(b)[2]; }
#define mult3(a,b,c) { (c)[0]=(a)*(b)[0]; (c)[1]=(a)*(b)[1]; (c)[2]=(a)*(b)[2]; }
#define multinv3(a,b,c) { (c)[0]=(b)[0]/(a); (c)[1]=(b)[1]/(a); (c)[2]=(b)[2]/(a); }
#define mult_plus3(a,b,c) { (c)[0]+=(a)*(b)[0]; (c)[1]+=(a)*(b)[1]; (c)[2]+=(a)*(b)[2]; }
#define copy3(a,c) { (c)[0]=(a)[0]; (c)[1]=(a)[1]; (c)[2]=(a)[2]; }
#define quaternion_product(a,b,c) { (c)[0]=a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3]; (c)[1]=a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2]; (c)[2]=a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1]; (c)[3]=a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0]; }
#define scale3(a,b,c) { (c)[0]=a*(c)[0]+(1-a)*(b)[0]; (c)[1]=a*(c)[1]+(1-a)*(b)[1]; (c)[2]=a*(c)[2]+(1-a)*(b)[2]; }
// calculate great-circle distance  on the sphere between points a and b given by decard coordinates.
// a and b assumed to have unit length.
realp dist_sphere_sq3(const real* a, const real* b);
#define for3(j) for(int j=0;j<3;j++)
#define random_real() ((real)rand()/(real)(RAND_MAX))
#define add_random_vector3(alpha, a, b) { (b)[0]=(a)[0]+(alpha)*(random_real()-0.5); (b)[1]=(a)[1]+(alpha)*(random_real()-0.5); (b)[2]=(a)[2]+(alpha)*(random_real()-0.5);  }

void add_random_cone3(real alpha, const real* A, real* B);
void matrixmult3(const real* a, const real* b, real* prod);
void invertmatrix3(const real* mat, real* inv);
realp distsq(int n, const real* a, const real* b);
realp dist_sphere_sq(int n, const real* a, const real* b);
realp normsq(int n, const real* a);
realp dot(int n, const real* a, const real* b);
void mult_sub(int n, real a, const real* b, real* c);
void mult_sub_ext(int n, real a, const real* b, const real* c, real* d);
void mult_add(int n, real a, const real* b, real* c);
void add_mult(int n, const real* a, real b, real* c); 
void sub_const(int n, real a, real* c);
void random_vector(int n, real* a);
void add_random_vector(real alpha, int n, const real* a, real* b);
void add_random_cone(real alpha, int n, const real* a, real* b);
void zero_vector(int n, real* a);
void copy_vector(int n, const real* a, real* b);
void sub_inplace(int n, const real* __restrict__ a, real* __restrict__ c);
void sub(int n, const real* __restrict__ a, const real* __restrict__ b, real* __restrict__ c);
void add_inplace(int n, const real* __restrict__ a, real* __restrict__ c);
void add_constant_inplace(int n, real a, real* __restrict__ c);
void negate_inplace(int n, real* __restrict__ c);
void negate_div(int n, real* __restrict__ a, real* __restrict__ c);
void const_div_inplace(int n, real a, real* __restrict__ c);
void linear_comb(int n, real a, const real* __restrict__ b, real c, real* __restrict__ d, real* __restrict__ e);

void gram_schmidt(int N, int K, real* __restrict__* V);

#endif