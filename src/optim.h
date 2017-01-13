#include "vector.h"

#define SDM_CONSTANT 0
#define SDM_INERTIAL 1
#define SDM_PROGR 2
// Compute minimum doing gradient descend.
// Minimizing function f(x) is of special form:
//   f(x)=<x|Hx>/2+<x|B>.
// Hessian part H is given by operator Q:x->H(x) with signature
//   void Q(int n, const real* x, real* y)
// Linear part is given as an operator L:y->y+B
//   void L(int n, real* y)
// The minimum is obtained as the limit of the sequence
// x(k+1)=x(k)-alpha(k)*grad_f(x(k))
// where grad_f(x) is the gradient of the minimizing function f(x):
//   grad_f(x)=Hx+B.
// Sequence alpha(k) is determined by mode:
//   SDM_CONSTANT : alpha(k)=mode_param
//   SDM_INERTIAL : alpha(k)=alpha(k-1)+mode_param*|grad_f(x(k-1))| if f(x(k))<f(x(k-1))
//                  alpha(k)=mode_param*|grad_f(x(k-1))| otherwise
// Initial approxiamtion is x(0)=a.
// Result is stored in the bufer 'a'.
// Length of the vector 'a' is 'n'.
// Progress can be displayed by function display:
//   void display(int d, int n, real* a, real* grad_f, real f, real res)
// where 
//   iter : iteration counter,
//   n : vector length,
//   a : current aproximation x(k),
//   grad_f : gradient of f(x) at a=x(k),
//   f : value of f(x(k)),
//   res : residual = |grad_f|.
// Return 0 on success (residual is less then epsilon)
//        1 if maximum number of iteration max_iter is reached
int steepest_descend(
  int n, real* a, 
  void (*TF)(const real* x, real* y, realp* E),
  int mode, real mode_param, real epsilon, int max_iter,
  void (*display)(int iter, real* a, real* grad_f, realp f, real res, real constres, real alpha, realp last_f, real last_grad, real err),
  real (*P)(real* a),
  void (*A)(const real* a, real* grad_f, real* f)
);

int flow_descend(
  int n, real* a, 
  void (*TF)(const real* x, real* y, realp* E),
  int mode, real mode_param, real epsilon, int max_iter,
  void (*display)(int iter, real* a, real* grad_f, realp f, real res, real constres, real alpha),
  real (*P)(real* a),
  void (*I)(int N, void (*F)(const real* x, real* g, realp* E), real T, real* X, realp* E, int* iter)
);

int lagrange_conjugate_quad(
  int N, int M, real* x0, 
  void (*Q)(const real* x, real* y),
  void (*L)(real* y),
  int mode, real mode_param, real epsilon, int max_iter,
  void (*display)(int iter, real* a, real* grad_f, realp f, real res, real alpha),
  void (*C)(const real* x, real* r),
  void (*D)(const real* x, const real* u, real* r),
  void (*P)(const real* x, const real* y, real* r),
  real initial_mu
 );

int lagrange_conjugate(
  int N, int M, real* x0, 
  void (*Q)(const real* x, real* y),
  void (*L)(real* y),
  int mode, real mode_param, real epsilon, int max_iter,
  void (*display)(int iter, real* a, real* grad_f, realp f, real res, real alpha),
  void (*C)(const real* x, real* r),
  void (*D)(const real* x, const real* u, real* r),
  void (*P)(const real* x, const real* y, real* r),
  real initial_mu
 );