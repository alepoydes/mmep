#include <math.h>

// считает длину вектора
#define normsq3(x) ((x)[0]*(x)[0]+(x)[1]*(x)[1]+(x)[2]*(x)[2])
// считает скалярное произведение векторов
#define dot3(x,y) (x)[0]*(y)[0]+(x)[1]*(y)[1]+(x)[2]*(y)[2]
//#define normalize3(x) { float t=(3-normsq3(x))*0.5; (x)[0]*=t; (x)[1]*=t; (x)[2]*=t; }
//#define normalize3(x) { float t=(1+normsq3(x))/2; (x)[0]/=t; (x)[1]/=t; (x)[2]/=t; }
#define normalize3(x) { double t=sqrt(normsq3(x)); (x)[0]/=t; (x)[1]/=t; (x)[2]/=t; }
#define middle3(x,y,z) { (z)[0]=((x)[0]+(y)[0])/2; (z)[1]=((x)[1]+(y)[1])/2; (z)[2]=((x)[2]+(y)[2])/2; normalize3(z); }
#define sub3(x,y,z) { (z)[0]=(x)[0]-(y)[0]; (z)[1]=(x)[1]-(y)[1]; (z)[2]=(x)[2]-(y)[2]; }
// проецирует y на касательное подпространство к x
#define tangent3(x,y) {	float t=dot3(x,y); (y)[0]-=t*(x)[0]; (y)[1]-=t*(x)[1]; (y)[2]-=t*(x)[2]; }
#define cross_minus3(a,b,c) { (c)[0]-=(a)[1]*(b)[2]-(a)[2]*(b)[1]; (c)[1]-=(a)[2]*(b)[0]-(a)[0]*(b)[2]; (c)[2]-=(a)[0]*(b)[1]-(a)[1]*(b)[0]; }
#define cross_plus3(a,b,c) { (c)[0]+=(a)[1]*(b)[2]-(a)[2]*(b)[1]; (c)[1]+=(a)[2]*(b)[0]-(a)[0]*(b)[2]; (c)[2]+=(a)[0]*(b)[1]-(a)[1]*(b)[0]; }
#define mult_minus3(a,b,c) { (c)[0]-=(a)*(b)[0]; (c)[1]-=(a)*(b)[1]; (c)[2]-=(a)*(b)[2]; }
#define mult_plus3(a,b,c) { (c)[0]+=(a)*(b)[0]; (c)[1]+=(a)*(b)[1]; (c)[2]+=(a)*(b)[2]; }

float normsq(int n, const float* a);
float dot(int n, const float* a, const float* b);
void mult_sub(int n, float a, const float* b, float* c);
void mult_add(int n, float a, const float* b, float* c);
void add_mult(int n, const float* a, float b, float* c); 
void sub_const(int n, float a, float* c);
void random_vector(int n, float* a);
void zero_vector(int n, float* a);
void copy_vector(int n, const float* a, float* b);
void sub_inplace(int n, const float* restrict a, float* restrict c);
void negate_inplace(int n, float* restrict c);